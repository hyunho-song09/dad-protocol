"""
dad.core.interaction — Stages 9-10: Interaction profiling, aggregation, and ranking.

Stage 9 — Interaction profiling:
    - 5 Å contact residue extraction (Bio.PDB, ported from GININA_Template.ipynb cell 19).
    - PLIP-based H-bond / hydrophobic / π-stacking classification (optional).

Stage 10 — Aggregation & ranking:
    - Long-format master CSV: proteins × ligands × pockets.
    - Composite rank score: z(vina)*w1 + z(cnn_pose)*w2 + z(cnn_affinity)*w3.
    - Heatmap data export.

Key defaults (from gnina.ipynb seed notebook)
----------------------------------------------
search_distance = 5.0 Å   (contact residue cutoff)

Implementation owner: Mr_Dock (Phase B).
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from dad.io import (
    DockingResult,
    InteractionProfile,
    ContactResidue,
    RankedPair,
)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 9 — Contact residue extraction
# ─────────────────────────────────────────────────────────────────────────────

def extract_contact_residues(
    receptor_pdb: str,
    docked_sdf: str,
    pose_index: int = 0,
    contact_distance: float = 5.0,
) -> List[ContactResidue]:
    """Extract protein residues within contact_distance of a docked ligand pose.

    Ports the Bio.PDB contact loop from GININA_Template.ipynb cell 19:
    - Parse receptor PDB with Bio.PDB.PDBParser.
    - Parse ligand coords from SDF (first model = best pose).
    - For each residue, find minimum atom-atom distance to any ligand atom.
    - Return residues with min_dist <= contact_distance, sorted by distance.
    """
    if not Path(receptor_pdb).exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")
    if not Path(docked_sdf).exists():
        raise FileNotFoundError(f"Docked SDF not found: {docked_sdf}")
    if contact_distance <= 0:
        return []

    ligand_coords = get_sdf_atom_coords(docked_sdf, model_index=pose_index)
    if len(ligand_coords) == 0:
        raise ValueError(f"No atom coordinates found in SDF model {pose_index}: {docked_sdf}")

    try:
        from Bio.PDB import PDBParser
    except ModuleNotFoundError:
        return _extract_contact_residues_plain(
            receptor_pdb=receptor_pdb,
            ligand_coords=ligand_coords,
            contact_distance=contact_distance,
        )

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", receptor_pdb)

    contacts: List[ContactResidue] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                min_dist = float("inf")
                for prot_atom in residue:
                    p_coord = prot_atom.get_coord()
                    dists = np.linalg.norm(ligand_coords - p_coord, axis=1)
                    d = float(dists.min())
                    if d < min_dist:
                        min_dist = d

                if min_dist <= contact_distance:
                    contacts.append(ContactResidue(
                        chain=chain.id,
                        res_id=int(residue.id[1]),
                        res_name=residue.resname.strip(),
                        min_dist=round(min_dist, 3),
                        interaction_type="contact",
                    ))

    contacts.sort(key=lambda r: r.min_dist)
    return contacts


def _extract_contact_residues_plain(
    receptor_pdb: str,
    ligand_coords: np.ndarray,
    contact_distance: float,
) -> List[ContactResidue]:
    """Fallback PDB parser used when BioPython is unavailable."""
    residue_min: Dict[Tuple[str, int, str], float] = {}
    for atom in _iter_pdb_atoms(receptor_pdb):
        p_coord = np.array([atom["x"], atom["y"], atom["z"]], dtype=float)
        d = float(np.linalg.norm(ligand_coords - p_coord, axis=1).min())
        key = (atom["chain"], atom["res_id"], atom["res_name"])
        if d <= contact_distance and (key not in residue_min or d < residue_min[key]):
            residue_min[key] = d

    contacts = [
        ContactResidue(
            chain=chain,
            res_id=res_id,
            res_name=res_name,
            min_dist=round(dist, 3),
            interaction_type="contact",
        )
        for (chain, res_id, res_name), dist in residue_min.items()
    ]
    contacts.sort(key=lambda r: r.min_dist)
    return contacts


def _iter_pdb_atoms(pdb_path: str) -> List[Dict[str, object]]:
    atoms: List[Dict[str, object]] = []
    for line in Path(pdb_path).read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            atoms.append({
                "chain": line[21].strip() or "A",
                "res_id": int(line[22:26]),
                "res_name": line[17:20].strip() or "UNK",
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })
        except ValueError:
            continue
    return atoms


def get_sdf_atom_coords(
    sdf_path: str,
    model_index: int = 0,
) -> np.ndarray:
    """Extract 3D atom coordinates from one model in a multi-model SDF file.

    Parses V2000 MOL block (same logic as GININA_Template.ipynb get_sdf_coords()).
    """
    if not Path(sdf_path).exists():
        raise FileNotFoundError(f"SDF not found: {sdf_path}")

    with open(sdf_path, "r") as fh:
        content = fh.read()

    models = content.split("$$$$")
    if model_index >= len(models):
        raise IndexError(
            f"model_index {model_index} out of range ({len(models)} models) in {sdf_path}"
        )

    lines = models[model_index].splitlines()
    coords = []
    start = None
    num_atoms = 0
    for i, line in enumerate(lines):
        if "V2000" in line:
            try:
                num_atoms = int(line.split()[0])
            except (ValueError, IndexError):
                num_atoms = 0
            start = i + 1
            break

    if start is None or num_atoms == 0:
        return np.zeros((0, 3))

    for i in range(start, start + num_atoms):
        if i >= len(lines):
            break
        parts = lines[i].split()
        if len(parts) >= 3:
            try:
                coords.append([float(parts[0]), float(parts[1]), float(parts[2])])
            except ValueError:
                pass

    return np.array(coords, dtype=float)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 9 — PLIP interaction classification
# ─────────────────────────────────────────────────────────────────────────────

def run_plip(
    receptor_pdb: str,
    docked_sdf: str,
    output_dir: str,
    plip_path: str = "plip",
    pose_index: int = 0,
) -> Optional[str]:
    """Run PLIP on a docked complex; return XML path or None if PLIP unavailable.

    Creates a merged receptor+ligand PDB (as in GININA_Template.ipynb cell 23),
    then calls PLIP to classify H-bonds, hydrophobic contacts, pi-stacking, etc.
    """
    try:
        # check plip is available
        subprocess.run([plip_path, "--help"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None

    # merge to a temp complex PDB
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    stem = Path(receptor_pdb).stem
    complex_pdb = str(out_path / f"{stem}_complex.pdb")

    _merge_pdb_sdf(receptor_pdb, docked_sdf, complex_pdb, pose_index=pose_index)

    xml_output = str(out_path / "report.xml")
    cmd = [plip_path, "-f", complex_pdb, "-x", "-o", str(out_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"PLIP failed: {proc.stderr[:1000]}")

    # PLIP writes report.xml in output_dir
    xml_candidates = list(out_path.glob("*.xml"))
    return str(xml_candidates[0]) if xml_candidates else None


def _merge_pdb_sdf(
    receptor_pdb: str,
    docked_sdf: str,
    output_pdb: str,
    pose_index: int = 0,
) -> None:
    """Merge receptor PDB + docked SDF pose into a single PDB (PLIP input)."""
    from rdkit import Chem

    with open(receptor_pdb, "r") as fh:
        rec_lines = [l for l in fh.readlines() if not l.startswith(("END", "CONECT"))]

    with open(docked_sdf, "r") as fh:
        content = fh.read()

    models = content.split("$$$$")
    if pose_index >= len(models):
        pose_index = 0
    mol_block = models[pose_index].strip() + "\n$$$$\n"

    suppl = Chem.SDMolSupplier()
    suppl.SetData(mol_block)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(f"Cannot parse SDF model {pose_index} from {docked_sdf}")

    lig_pdb_block = Chem.MolToPDBBlock(mol)
    lig_lines = [
        l + "\n" for l in lig_pdb_block.splitlines()
        if l.startswith(("ATOM", "HETATM", "CONECT"))
    ]

    with open(output_pdb, "w") as fh:
        fh.writelines(rec_lines)
        fh.writelines(lig_lines)
        fh.write("END\n")


def profile_interactions(
    docking_results: List[DockingResult],
    structures: Dict[str, object],
    config: Dict,
    output_dir: str = "results/interactions",
) -> List[InteractionProfile]:
    """Generate full interaction profiles for all docking results."""
    int_cfg = config.get("interaction", {})
    contact_distance = int_cfg.get("contact_distance", 5.0)
    plip_path = int_cfg.get("plip_path", "plip")
    use_plip = int_cfg.get("use_plip", True)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    profiles: List[InteractionProfile] = []
    for result in docking_results:
        struct = structures.get(result.seq_id)
        if struct is None:
            continue

        receptor_pdb = getattr(struct, "pdb_path", None) or str(struct)
        pair_dir = str(out_path / f"{result.seq_id}__{result.lig_id}__p{result.pocket_rank}")

        contacts = extract_contact_residues(
            receptor_pdb=receptor_pdb,
            docked_sdf=result.output_sdf,
            pose_index=0,
            contact_distance=contact_distance,
        )

        plip_xml = None
        n_hbonds = None
        n_hydrophobic = None
        if use_plip:
            try:
                plip_xml = run_plip(
                    receptor_pdb=receptor_pdb,
                    docked_sdf=result.output_sdf,
                    output_dir=pair_dir,
                    plip_path=plip_path,
                )
                if plip_xml:
                    n_hbonds, n_hydrophobic = _parse_plip_xml(plip_xml)
            except Exception:
                pass

        profiles.append(InteractionProfile(
            seq_id=result.seq_id,
            lig_id=result.lig_id,
            pocket_rank=result.pocket_rank,
            pose_rank=result.best_pose.pose_rank,
            contact_residues=contacts,
            n_hbonds=n_hbonds,
            n_hydrophobic=n_hydrophobic,
            plip_xml_path=plip_xml,
        ))

    return profiles


def _parse_plip_xml(xml_path: str) -> Tuple[Optional[int], Optional[int]]:
    """Extract hbond and hydrophobic counts from PLIP XML."""
    try:
        import xml.etree.ElementTree as ET
        tree = ET.parse(xml_path)
        root = tree.getroot()
        hbonds = len(root.findall(".//hbond"))
        hydros = len(root.findall(".//hydrophobic_interaction"))
        return hbonds, hydros
    except Exception:
        return None, None


# ─────────────────────────────────────────────────────────────────────────────
# Stage 10 — Aggregation & ranking
# ─────────────────────────────────────────────────────────────────────────────

def aggregate_results(
    docking_results: List[DockingResult],
    interaction_profiles: List[InteractionProfile],
    config: Dict,
) -> List[RankedPair]:
    """Aggregate docking and interaction results into a ranked master table.

    Composite ranking:
        z_vina      = z-score of (-vina_affinity)   (negate: lower vina = better)
        z_cnn_pose  = z-score of cnn_pose_score      (higher = better)
        z_cnn_aff   = z-score of cnn_affinity        (higher = better)
        composite   = w_vina * z_vina + w_cnn_pose * z_cnn_pose + w_cnn_aff * z_cnn_aff
    """
    rank_cfg = config.get("ranking", {})
    weights = rank_cfg.get("score_weights", {"vina": 1.0, "cnn_pose": 1.0, "cnn_affinity": 1.0})
    w_vina = weights.get("vina", 1.0)
    w_pose = weights.get("cnn_pose", 1.0)
    w_aff = weights.get("cnn_affinity", 1.0)

    profile_map: Dict[Tuple, InteractionProfile] = {
        (p.seq_id, p.lig_id, p.pocket_rank): p
        for p in interaction_profiles
    }

    vina_vals = np.array([-r.best_pose.vina_affinity for r in docking_results])
    pose_vals = np.array([r.best_pose.cnn_pose_score for r in docking_results])
    aff_vals = np.array([r.best_pose.cnn_affinity for r in docking_results])

    z_vina = _zscore(vina_vals)
    z_pose = _zscore(pose_vals)
    z_aff = _zscore(aff_vals)
    composite = w_vina * z_vina + w_pose * z_pose + w_aff * z_aff

    rows = []
    for i, result in enumerate(docking_results):
        profile = profile_map.get((result.seq_id, result.lig_id, result.pocket_rank))
        n_contacts = len(profile.contact_residues) if profile else None
        n_hbonds = profile.n_hbonds if profile else None

        rows.append({
            "seq_id": result.seq_id,
            "lig_id": result.lig_id,
            "pocket_rank": result.pocket_rank,
            "vina_affinity": result.best_pose.vina_affinity,
            "cnn_pose_score": result.best_pose.cnn_pose_score,
            "cnn_affinity": result.best_pose.cnn_affinity,
            "composite_rank_score": float(composite[i]),
            "n_contacts": n_contacts,
            "n_hbonds": n_hbonds,
        })

    rows.sort(key=lambda r: r["composite_rank_score"], reverse=True)
    ranked = []
    for rank, row in enumerate(rows, start=1):
        ranked.append(RankedPair(
            seq_id=row["seq_id"],
            lig_id=row["lig_id"],
            pocket_rank=row["pocket_rank"],
            vina_affinity=row["vina_affinity"],
            cnn_pose_score=row["cnn_pose_score"],
            cnn_affinity=row["cnn_affinity"],
            composite_rank_score=row["composite_rank_score"],
            overall_rank=rank,
            n_contacts=row["n_contacts"],
            n_hbonds=row["n_hbonds"],
        ))
    return ranked


def _zscore(arr: np.ndarray) -> np.ndarray:
    """Z-score normalise; returns zeros if std == 0."""
    if len(arr) == 0:
        return arr
    std = arr.std()
    if std == 0:
        return np.zeros_like(arr)
    return (arr - arr.mean()) / std


def write_master_csv(
    ranked_pairs: List[RankedPair],
    output_csv: str,
) -> None:
    """Write the ranked master table to a long-format CSV file.

    Columns match the Revision.txt Table X format:
    seq_id, lig_id, pocket_rank, vina_affinity, cnn_pose_score,
    cnn_affinity, composite_rank_score, overall_rank, n_contacts, n_hbonds
    """
    Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "overall_rank", "seq_id", "lig_id", "pocket_rank",
        "vina_affinity", "cnn_pose_score", "cnn_affinity",
        "composite_rank_score", "n_contacts", "n_hbonds",
    ]
    with open(output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for rp in ranked_pairs:
            writer.writerow({
                "overall_rank": rp.overall_rank,
                "seq_id": rp.seq_id,
                "lig_id": rp.lig_id,
                "pocket_rank": rp.pocket_rank,
                "vina_affinity": rp.vina_affinity,
                "cnn_pose_score": rp.cnn_pose_score,
                "cnn_affinity": rp.cnn_affinity,
                "composite_rank_score": rp.composite_rank_score,
                "n_contacts": rp.n_contacts,
                "n_hbonds": rp.n_hbonds,
            })


def build_heatmap_matrix(
    ranked_pairs: List[RankedPair],
    score_column: str = "cnn_affinity",
) -> Tuple[List[str], List[str], np.ndarray]:
    """Build a 2D matrix for heatmap visualization.

    Rows = proteins (seq_id), columns = ligands (lig_id).
    Values = best score across all pockets for each (protein, ligand) pair.
    """
    valid_cols = {"cnn_affinity", "vina_affinity", "composite_rank_score"}
    if score_column not in valid_cols:
        raise ValueError(f"score_column must be one of {valid_cols}, got {score_column!r}")

    protein_ids = sorted({rp.seq_id for rp in ranked_pairs})
    ligand_ids = sorted({rp.lig_id for rp in ranked_pairs})
    p_idx = {p: i for i, p in enumerate(protein_ids)}
    l_idx = {l: i for i, l in enumerate(ligand_ids)}

    # initialise with NaN; fill with best-pocket score per pair
    matrix = np.full((len(protein_ids), len(ligand_ids)), np.nan)
    for rp in ranked_pairs:
        val = getattr(rp, score_column)
        pi = p_idx[rp.seq_id]
        li = l_idx[rp.lig_id]
        if np.isnan(matrix[pi, li]) or val > matrix[pi, li]:
            matrix[pi, li] = val

    return protein_ids, ligand_ids, matrix
