"""
dad.core.docking — Stages 6-8: Ligand preparation, box config, and GNINA docking.

Stage 6 — Ligand preparation: RDKit 3D embedding + MMFF94 minimization,
           Open Babel protonation at pH 7.4, optional tautomer/stereo enumeration.
Stage 7 — Auto box configuration: pocket centroid + DAD box-sizing rule.
Stage 8 — GNINA 1.3.2 docking: all (protein × ligand × pocket) combinations.

Key defaults (from gnina.ipynb seed notebook)
----------------------------------------------
exhaustiveness = 32
num_modes      = 9
seed           = 0
size_x/y/z     = 22 (minimum; auto-scaled to ligand dims)
contact_cutoff = 5.0 Å  (used in Stage 9, defined here for reference)

Implementation owner: Mr_Dock (Phase B).
"""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from dad.io import (
    ProteinInput,
    StructurePrediction,
    PocketResult,
    LigandInput,
    PreparedLigand,
    DockingBox,
    DockingResult,
    DockingPose,
)

# ─────────────────────────────────────────────────────────────────────────────
# Stage 6 — Ligand preparation
# ─────────────────────────────────────────────────────────────────────────────

def prepare_ligands(
    ligands: List[LigandInput],
    config: Dict,
    output_dir: str = "results/ligands",
) -> List[PreparedLigand]:
    """Prepare all input SMILES as 3D docking-ready SDF files."""
    lig_cfg = config.get("ligand", {})
    ph = lig_cfg.get("ph", 7.4)
    force_field = lig_cfg.get("force_field", "MMFF94")
    enumerate_tautomers = lig_cfg.get("enumerate_tautomers", False)
    enumerate_stereo = lig_cfg.get("enumerate_stereo", False)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    results = []
    for lig in ligands:
        sdf_path = str(out_path / f"{lig.lig_id}.sdf")
        try:
            prepared = smiles_to_3d_sdf(
                smiles=lig.smiles,
                lig_id=lig.lig_id,
                output_sdf=sdf_path,
                force_field=force_field,
                ph=ph,
                enumerate_tautomers=enumerate_tautomers,
                enumerate_stereo=enumerate_stereo,
            )
            results.append(prepared)
        except Exception as exc:
            import logging
            logging.getLogger(__name__).warning("prepare_ligands: %s failed — %s", lig.lig_id, exc)
    return results


def smiles_to_3d_sdf(
    smiles: str,
    lig_id: str,
    output_sdf: str,
    force_field: str = "MMFF94",
    ph: float = 7.4,
    add_hydrogens: bool = True,
    enumerate_tautomers: bool = False,
    enumerate_stereo: bool = False,
    random_seed: int = 0,
) -> PreparedLigand:
    """Convert a SMILES string to a 3D-optimised SDF file.

    Uses Open Babel (obabel) for protonation and 3D generation to match
    the seed notebook workflow (GININA_Template.ipynb cell 9-10).
    Falls back to RDKit embedding if obabel is unavailable.
    """
    try:
        return _smiles_to_3d_obabel(smiles, lig_id, output_sdf, ph, random_seed)
    except (FileNotFoundError, RuntimeError):
        return _smiles_to_3d_rdkit(
            smiles, lig_id, output_sdf, force_field,
            add_hydrogens, enumerate_tautomers, enumerate_stereo, random_seed
        )


def _smiles_to_3d_obabel(
    smiles: str,
    lig_id: str,
    output_sdf: str,
    ph: float = 7.4,
    random_seed: int = 0,
) -> PreparedLigand:
    """3D generation via obabel (mirrors GININA_Template.ipynb cells 9-10)."""
    cmd = [
        "obabel", f"-:{smiles}",
        "-O", output_sdf,
        "--gen3d",
        "-p", str(ph),
        "--partialcharge", "gasteiger",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not Path(output_sdf).exists():
        raise RuntimeError(f"obabel failed: {result.stderr.strip()}")

    max_dim = compute_ligand_max_dim(output_sdf)
    mol_weight = _get_mol_weight_obabel(output_sdf)

    return PreparedLigand(
        lig_id=lig_id,
        sdf_path=output_sdf,
        smiles_canonical=smiles,
        mol_weight=mol_weight,
        max_dim=max_dim,
    )


def _smiles_to_3d_rdkit(
    smiles: str,
    lig_id: str,
    output_sdf: str,
    force_field: str = "MMFF94",
    add_hydrogens: bool = True,
    enumerate_tautomers: bool = False,
    enumerate_stereo: bool = False,
    random_seed: int = 0,
) -> PreparedLigand:
    """3D generation via RDKit (fallback path)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse SMILES: {smiles!r}")

    canonical_smiles = Chem.MolToSmiles(mol)
    mol_weight = Descriptors.MolWt(mol)

    if add_hydrogens:
        mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    ret = AllChem.EmbedMolecule(mol, params)
    if ret == -1:
        # try ETDG if ETKDGv3 fails
        ret = AllChem.EmbedMolecule(mol, AllChem.ETDG())
    if ret == -1:
        raise RuntimeError(f"3D embedding failed for SMILES: {smiles!r}")

    if force_field.upper() == "MMFF94":
        ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    else:
        ff = AllChem.UFFGetMoleculeForceField(mol)
    if ff is not None:
        ff.Minimize()

    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()

    max_dim = compute_ligand_max_dim(output_sdf)
    return PreparedLigand(
        lig_id=lig_id,
        sdf_path=output_sdf,
        smiles_canonical=canonical_smiles,
        mol_weight=mol_weight,
        max_dim=max_dim,
    )


def _get_mol_weight_obabel(sdf_path: str) -> float:
    """Read molecular weight from SDF property block, else use RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol:
            return Descriptors.MolWt(mol)
    except Exception:
        pass
    return 0.0


def compute_ligand_max_dim(sdf_path: str) -> float:
    """Compute the maximum pairwise atom distance in a 3D SDF conformer.

    Used by Stage 7 for auto box sizing:
        box_side = max(22 Å, max_dim + 10 Å)

    Mirrors the approach in GININA_Template.ipynb get_sdf_coords() logic.
    """
    import numpy as np

    coords = _parse_sdf_v2000_coords(sdf_path, model_index=0)
    if len(coords) < 2:
        return 0.0

    # efficient pairwise max via broadcasting on small molecules
    diff = coords[:, None, :] - coords[None, :, :]
    dists = np.sqrt((diff ** 2).sum(axis=-1))
    return float(dists.max())


def _parse_sdf_v2000_coords(sdf_path: str, model_index: int = 0) -> "np.ndarray":
    """Parse atom xyz from V2000 SDF (same logic as GININA_Template.ipynb cell 19)."""
    import numpy as np

    with open(sdf_path, "r") as fh:
        content = fh.read()

    models = content.split("$$$$")
    if model_index >= len(models):
        raise IndexError(f"model_index {model_index} out of range ({len(models)} models)")

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
# Stage 7 — Auto box configuration
# ─────────────────────────────────────────────────────────────────────────────

def build_docking_boxes(
    pockets: Dict[str, List[PocketResult]],
    ligands: List[PreparedLigand],
    config: Dict,
) -> List[DockingBox]:
    """Build GNINA docking boxes for all (protein × ligand × pocket) triples.

    Box sizing rule (zero-tuning, DAD Stage 7):
        side = max(config.box_size_min, ligand.max_dim + config.box_padding)
    """
    dock_cfg = config.get("docking", {})
    box_size_min = dock_cfg.get("box_size_min", 22.0)
    box_padding = dock_cfg.get("box_padding", 10.0)

    boxes = []
    for lig in ligands:
        side = max(box_size_min, lig.max_dim + box_padding)
        for seq_id, pocket_list in pockets.items():
            for pocket in pocket_list:
                boxes.append(DockingBox(
                    seq_id=seq_id,
                    lig_id=lig.lig_id,
                    pocket_rank=pocket.pocket_rank,
                    center_x=pocket.center_x,
                    center_y=pocket.center_y,
                    center_z=pocket.center_z,
                    size_x=side,
                    size_y=side,
                    size_z=side,
                ))
    return boxes


# ─────────────────────────────────────────────────────────────────────────────
# Stage 8 — GNINA docking
# ─────────────────────────────────────────────────────────────────────────────

def run_docking_batch(
    structures: Dict[str, StructurePrediction],
    ligands: Dict[str, PreparedLigand],
    boxes: List[DockingBox],
    config: Dict,
    output_dir: str = "results/docking",
) -> List[DockingResult]:
    """Run GNINA docking for all (protein × ligand × pocket) combinations."""
    dock_cfg = config.get("docking", {})
    gnina_path = dock_cfg.get("gnina_path", "gnina")
    exhaustiveness = dock_cfg.get("exhaustiveness", 32)
    num_modes = dock_cfg.get("num_modes", 9)
    seed = dock_cfg.get("seed", 0)
    cnn_scoring = dock_cfg.get("cnn_scoring", "rescore")
    addH = dock_cfg.get("addH", True)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    results = []
    for box in boxes:
        structure = structures.get(box.seq_id)
        ligand = ligands.get(box.lig_id)
        if structure is None or ligand is None:
            continue

        out_sdf = str(out_path / f"{box.seq_id}__{box.lig_id}__pocket{box.pocket_rank}.sdf")
        try:
            result = run_gnina_single(
                receptor_pdb=structure.pdb_path,
                ligand_sdf=ligand.sdf_path,
                box=box,
                output_sdf=out_sdf,
                gnina_path=gnina_path,
                exhaustiveness=exhaustiveness,
                num_modes=num_modes,
                seed=seed,
                cnn_scoring=cnn_scoring,
                addH=addH,
            )
            results.append(result)
        except Exception as exc:
            import logging
            logging.getLogger(__name__).warning(
                "run_docking_batch: %s × %s pocket%d failed — %s",
                box.seq_id, box.lig_id, box.pocket_rank, exc
            )
    return results


def run_gnina_single(
    receptor_pdb: str,
    ligand_sdf: str,
    box: DockingBox,
    output_sdf: str,
    gnina_path: str = "gnina",
    exhaustiveness: int = 32,
    num_modes: int = 9,
    seed: int = 0,
    cnn_scoring: str = "rescore",
    addH: bool = True,
) -> DockingResult:
    """Run GNINA for a single (receptor, ligand, box) combination.

    Command mirrors GININA_Template.ipynb cell 12:
        gnina -r structure.pdb -l ala_ile.sdf
              --center_x -18.1706 --center_y 3.1775 --center_z 7.1226
              --size_x 22 --size_y 22 --size_z 22
              --exhaustiveness 32 --num_modes 9 --seed 0
              --out docked_ala_ile.sdf
    """
    if not Path(receptor_pdb).exists():
        raise FileNotFoundError(f"Receptor not found: {receptor_pdb}")
    if not Path(ligand_sdf).exists():
        raise FileNotFoundError(f"Ligand SDF not found: {ligand_sdf}")

    cmd = [
        gnina_path,
        "-r", receptor_pdb,
        "-l", ligand_sdf,
        "--center_x", str(box.center_x),
        "--center_y", str(box.center_y),
        "--center_z", str(box.center_z),
        "--size_x", str(box.size_x),
        "--size_y", str(box.size_y),
        "--size_z", str(box.size_z),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--seed", str(seed),
        "--cnn_scoring", cnn_scoring,
        "--out", output_sdf,
    ]
    if addH:
        cmd.append("--addH")

    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0

    if proc.returncode != 0:
        raise RuntimeError(
            f"GNINA failed (exit {proc.returncode}):\n{proc.stderr[:2000]}"
        )

    result = parse_gnina_output_sdf(
        output_sdf=output_sdf,
        seq_id=box.seq_id,
        lig_id=box.lig_id,
        pocket_rank=box.pocket_rank,
    )
    result.elapsed_seconds = elapsed

    # capture version from stderr banner
    ver_match = re.search(r"gnina\s+v(\S+)", proc.stderr + proc.stdout)
    if ver_match:
        result.gnina_version = ver_match.group(1)

    return result


def parse_gnina_output_sdf(
    output_sdf: str,
    seq_id: str,
    lig_id: str,
    pocket_rank: int,
) -> DockingResult:
    """Parse GNINA multi-model SDF output into a DockingResult.

    GNINA writes scores as REMARK lines per model:
        REMARK minimizedAffinity <vina_affinity>
        REMARK CNNscore <cnn_pose_score>
        REMARK CNNaffinity <cnn_affinity>
    """
    if not Path(output_sdf).exists():
        raise FileNotFoundError(f"GNINA output SDF not found: {output_sdf}")

    with open(output_sdf, "r") as fh:
        content = fh.read()

    models = [m for m in content.split("$$$$") if m.strip()]
    if not models:
        raise ValueError(f"No models in GNINA output: {output_sdf}")

    poses: List[DockingPose] = []
    for rank_0, model_text in enumerate(models):
        vina = _parse_remark_float(model_text, "minimizedAffinity")
        cnn_pose = _parse_remark_float(model_text, "CNNscore")
        cnn_aff = _parse_remark_float(model_text, "CNNaffinity")

        # fallback: try VINA RESULT line
        if vina is None:
            m = re.search(r"REMARK VINA RESULT:\s+([-\d.]+)", model_text)
            if m:
                vina = float(m.group(1))

        poses.append(DockingPose(
            pose_rank=rank_0 + 1,
            vina_affinity=vina if vina is not None else 0.0,
            cnn_pose_score=cnn_pose if cnn_pose is not None else 0.0,
            cnn_affinity=cnn_aff if cnn_aff is not None else 0.0,
        ))

    best = poses[0]
    return DockingResult(
        seq_id=seq_id,
        lig_id=lig_id,
        pocket_rank=pocket_rank,
        output_sdf=output_sdf,
        poses=poses,
        best_pose=best,
    )


def _parse_remark_float(model_text: str, keyword: str) -> Optional[float]:
    pattern = rf"REMARK\s+{re.escape(keyword)}\s+([-\d.]+)"
    m = re.search(pattern, model_text)
    return float(m.group(1)) if m else None


# ─────────────────────────────────────────────────────────────────────────────
# Tier 1 replay mode — no GPU required
# ─────────────────────────────────────────────────────────────────────────────

# Hardcoded ground truth from Revision.txt Table X.
# Key: (target_label, ligand_id) → (vina_affinity, cnn_pose_score, cnn_affinity)
# target_label matches the TARGET_ALIASES keys (case-insensitive prefix match).
_REVISION_GROUND_TRUTH: Dict[Tuple[str, str], Tuple[float, float, float]] = {
    ("MCP", "AlaIle"):  (-5.61, 0.8995, 4.416),
    ("MCP", "GlyVal"):  (-5.18, 0.8468, 4.022),
    ("Crp", "AlaIle"):  (-6.01, 0.6228, 4.642),
    ("Crp", "GlyVal"):  (-5.59, 0.7353, 4.574),
    ("RbsB", "AlaIle"): (-4.92, 0.4447, 3.516),
    ("RbsB", "GlyVal"): (-4.44, 0.9144, 3.896),
}

# seq_id / case_id normalisation table — maps variant names to canonical keys.
_TARGET_ALIASES: Dict[str, str] = {
    "MCP":               "MCP",
    "NA23_RS01195":      "MCP",
    "mcpB":              "MCP",
    "mcp":               "MCP",
    "MCP sensory domain":"MCP",
    "Crp":               "Crp",
    "CRP":               "Crp",
    "NA23_RS08105":      "Crp",
    "crp":               "Crp",
    "Crp/Fnr family regulator": "Crp",
    "RbsB":              "RbsB",
    "rbsB":              "RbsB",
    "NA23_RS00870":      "RbsB",
    "rbs":               "RbsB",
    "RbsB (substrate-binding protein)": "RbsB",
}

_LIGAND_ALIASES: Dict[str, str] = {
    "AlaIle":  "AlaIle",
    "Ala-Ile": "AlaIle",
    "ala_ile": "AlaIle",
    "ala-ile": "AlaIle",
    "CC[C@H](C)[C@@H](C(=O)O)NC(=O)[C@H](C)N": "AlaIle",
    "GlyVal":  "GlyVal",
    "Gly-Val": "GlyVal",
    "gly_val": "GlyVal",
    "gly-val": "GlyVal",
    "CC(C)[C@@H](C(=O)O)NC(=O)CN": "GlyVal",
}


def _normalise_target(case_id: str) -> Optional[str]:
    """Return canonical target key for a case_id/seq_id string, or None."""
    if case_id in _TARGET_ALIASES:
        return _TARGET_ALIASES[case_id]
    # prefix match for long names
    lower = case_id.lower()
    for alias, canonical in _TARGET_ALIASES.items():
        if lower.startswith(alias.lower()) or alias.lower().startswith(lower):
            return canonical
    return None


def _normalise_ligand(ligand_id: str) -> Optional[str]:
    """Return canonical ligand key for a ligand_id/smiles string, or None."""
    if ligand_id in _LIGAND_ALIASES:
        return _LIGAND_ALIASES[ligand_id]
    lower = ligand_id.lower()
    for alias, canonical in _LIGAND_ALIASES.items():
        if lower == alias.lower():
            return canonical
    return None


def _parse_revision_txt(revision_table_path: Path) -> Dict[Tuple[str, str], Tuple[float, float, float]]:
    """Parse Revision.txt Table X and return ground truth dict.

    Accepts lines of the form (tab-separated):
        <Target>    <Ligand>    <Vina>    <CNN pose>    <CNN affinity>    [<Product>]

    Falls back to the hardcoded _REVISION_GROUND_TRUTH if file cannot be parsed.
    """
    table: Dict[Tuple[str, str], Tuple[float, float, float]] = {}
    try:
        text = revision_table_path.read_text(encoding="utf-8")
    except OSError:
        return {}

    for line in text.splitlines():
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 5:
            continue
        target_raw, lig_raw = parts[0], parts[1]
        try:
            vina = float(parts[2])
            cnn_pose = float(parts[3])
            cnn_aff = float(parts[4])
        except ValueError:
            continue
        t_key = _normalise_target(target_raw)
        l_key = _normalise_ligand(lig_raw)
        if t_key and l_key:
            table[(t_key, l_key)] = (vina, cnn_pose, cnn_aff)
    return table


def replay_from_ground_truth(
    case_id: str,
    revision_table_path: Path,
    ligand_id: str,
    pocket_rank: int = 1,
) -> DockingResult:
    """Construct DockingResult from precomputed Revision.txt ground truth.

    For Tier 1 replay mode (no GPU). Sets gnina_version='precomputed_revision'
    as the pose source marker so downstream stages can distinguish replay from
    live GNINA runs.

    Used when execution.mode == 'tier1_replay' to validate pipeline I/O shape
    without running GNINA.

    Parameters
    ----------
    case_id : str
        Target identifier — any variant of MCP/Crp/RbsB accepted
        (e.g. "NA23_RS01195", "MCP", "Crp/Fnr family regulator").
    revision_table_path : Path
        Path to Revision.txt containing Table X (tab-separated).
        If the file cannot be parsed, falls back to the hardcoded ground truth.
    ligand_id : str
        Ligand identifier — any variant of AlaIle/GlyVal accepted
        (e.g. "Ala-Ile", "GlyVal", "gly_val").
    pocket_rank : int, optional
        Pocket rank to record (default 1; Tier 1 always uses top pocket).

    Returns
    -------
    DockingResult
        Populated from ground truth scores. ``output_sdf`` is set to the
        sentinel string ``"<replay:no_sdf>"`` — callers in replay mode must
        not attempt to open this path.

    Raises
    ------
    KeyError
        If case_id / ligand_id cannot be resolved to a known Tier 1 pair.
    """
    t_key = _normalise_target(case_id)
    l_key = _normalise_ligand(ligand_id)

    if t_key is None:
        raise KeyError(
            f"replay_from_ground_truth: unknown target '{case_id}'. "
            f"Known targets: {sorted(set(_TARGET_ALIASES.values()))}"
        )
    if l_key is None:
        raise KeyError(
            f"replay_from_ground_truth: unknown ligand '{ligand_id}'. "
            f"Known ligands: {sorted(set(_LIGAND_ALIASES.values()))}"
        )

    # try to parse from file; fall back to hardcoded table
    file_table = _parse_revision_txt(revision_table_path)
    gt_table = file_table if file_table else _REVISION_GROUND_TRUTH

    gt_key = (t_key, l_key)
    if gt_key not in gt_table:
        # try hardcoded as final fallback regardless
        if gt_key not in _REVISION_GROUND_TRUTH:
            raise KeyError(
                f"replay_from_ground_truth: no ground truth for ({t_key!r}, {l_key!r}). "
                f"Available pairs: {sorted(_REVISION_GROUND_TRUTH.keys())}"
            )
        vina, cnn_pose, cnn_aff = _REVISION_GROUND_TRUTH[gt_key]
    else:
        vina, cnn_pose, cnn_aff = gt_table[gt_key]

    best = DockingPose(
        pose_rank=1,
        vina_affinity=vina,
        cnn_pose_score=cnn_pose,
        cnn_affinity=cnn_aff,
    )
    return DockingResult(
        seq_id=case_id,
        lig_id=ligand_id,
        pocket_rank=pocket_rank,
        output_sdf="<replay:no_sdf>",
        poses=[best],
        best_pose=best,
        # gnina_version carries the pose source marker for replay mode
        gnina_version="precomputed_revision",
        elapsed_seconds=0.0,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Reference-free pose selection (Blocker 4 / Nature Protocols Phase C)
# ─────────────────────────────────────────────────────────────────────────────

def select_consensus_pose(
    poses_sdf: Path,
    receptor_pdb: Path,
    method: str = "ensemble",
    rmsd_cluster_threshold: float = 2.0,
) -> dict:
    """Select a single pose from a multi-model GNINA SDF without using a crystal reference.

    Parameters
    ----------
    poses_sdf : Path
        Path to a multi-model GNINA SDF (e.g. docked.sdf with 9 models).
    receptor_pdb : Path
        Path to receptor PDB (used by PLIP method; otherwise read for validation only).
    method : str
        One of: ``ensemble``, ``cluster``, ``posebusters``, ``plip``, ``consensus``.
    rmsd_cluster_threshold : float
        Å threshold for RMSD-based clustering (used by ``cluster`` and ``consensus`` methods).

    Returns
    -------
    dict with keys:
        ``selected_pose_index`` — 0-based index into the model list.
        ``score_breakdown``     — dict with per-method subscores used for selection.
        ``method_used``         — method that actually made the selection
                                  (may differ from ``method`` when fallback occurs).
    """
    poses_sdf = Path(poses_sdf)
    if not poses_sdf.exists():
        raise FileNotFoundError(f"poses_sdf not found: {poses_sdf}")

    raw_poses = _read_gnina_poses(poses_sdf)
    if not raw_poses:
        raise ValueError(f"No poses parsed from {poses_sdf}")

    if method == "ensemble":
        return _select_ensemble(raw_poses)
    elif method == "cluster":
        return _select_cluster(raw_poses, poses_sdf, rmsd_cluster_threshold)
    elif method == "posebusters":
        return _select_posebusters(raw_poses, poses_sdf, receptor_pdb)
    elif method == "plip":
        return _select_plip(raw_poses, poses_sdf, receptor_pdb)
    elif method == "consensus":
        return _select_consensus(raw_poses, poses_sdf, receptor_pdb, rmsd_cluster_threshold)
    else:
        raise ValueError(
            f"Unknown method '{method}'. Choose from: ensemble, cluster, posebusters, plip, consensus."
        )


# ── internal helpers ──────────────────────────────────────────────────────────

def _read_gnina_poses(poses_sdf: Path) -> List[Dict]:
    """Parse all models from a GNINA SDF into a list of score dicts."""
    content = poses_sdf.read_text(encoding="utf-8", errors="replace")
    models = [m for m in content.split("$$$$") if m.strip()]
    result = []
    for i, model_text in enumerate(models):
        vina    = _parse_remark_float(model_text, "minimizedAffinity")
        cnn_pose = _parse_remark_float(model_text, "CNNscore")
        cnn_aff  = _parse_remark_float(model_text, "CNNaffinity")
        result.append({
            "index":        i,
            "model_text":   model_text,
            "vina":         vina    if vina    is not None else 0.0,
            "cnn_pose":     cnn_pose if cnn_pose is not None else 0.0,
            "cnn_affinity": cnn_aff  if cnn_aff  is not None else 0.0,
        })
    return result


def _ensemble_score(pose: Dict) -> float:
    """CNN_pose × 0.6 + CNN_affinity × 0.4 (higher is better)."""
    return pose["cnn_pose"] * 0.6 + pose["cnn_affinity"] * 0.4


def _select_ensemble(raw_poses: List[Dict]) -> dict:
    """Select by CNN ensemble score (reference-free, always available)."""
    scored = [(i, _ensemble_score(p)) for i, p in enumerate(raw_poses)]
    best_i, best_score = max(scored, key=lambda x: x[1])
    return {
        "selected_pose_index": best_i,
        "score_breakdown": {
            "method": "ensemble",
            "ensemble_scores": [s for _, s in scored],
            "best_ensemble_score": best_score,
        },
        "method_used": "ensemble",
    }


def _select_cluster(
    raw_poses: List[Dict],
    poses_sdf: Path,
    threshold: float = 2.0,
) -> dict:
    """Select the representative of the largest RMSD cluster.

    Uses RDKit GetBestRMS for pairwise ligand RMSD. Falls back to ensemble
    if RDKit is unavailable or fewer than 2 poses are present.
    """
    if len(raw_poses) < 2:
        result = _select_ensemble(raw_poses)
        result["method_used"] = "cluster→ensemble(single_pose)"
        return result

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        import warnings
        warnings.warn(
            "select_consensus_pose(method='cluster'): RDKit not available — falling back to ensemble.",
            stacklevel=3,
        )
        result = _select_ensemble(raw_poses)
        result["method_used"] = "cluster→ensemble(no_rdkit)"
        return result

    # Load all poses as RDKit molecules
    suppl = Chem.SDMolSupplier(str(poses_sdf), removeHs=False, sanitize=False)
    mols = [m for m in suppl if m is not None]
    if len(mols) < 2:
        result = _select_ensemble(raw_poses)
        result["method_used"] = "cluster→ensemble(rdkit_load_failed)"
        return result

    n = len(mols)
    # pairwise RMSD matrix
    import numpy as np
    rmsd_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            try:
                r = AllChem.GetBestRMS(mols[i], mols[j])
            except Exception:
                r = 999.0
            rmsd_matrix[i, j] = r
            rmsd_matrix[j, i] = r

    # greedy clustering: group poses within threshold
    assigned = [-1] * n
    cluster_id = 0
    for i in range(n):
        if assigned[i] >= 0:
            continue
        assigned[i] = cluster_id
        for j in range(i + 1, n):
            if assigned[j] < 0 and rmsd_matrix[i, j] <= threshold:
                assigned[j] = cluster_id
        cluster_id += 1

    # largest cluster
    from collections import Counter
    cluster_sizes = Counter(assigned)
    largest_cluster = cluster_sizes.most_common(1)[0][0]
    members = [i for i, c in enumerate(assigned) if c == largest_cluster]

    # within largest cluster pick highest ensemble score
    best_i = max(members, key=lambda i: _ensemble_score(raw_poses[i]))

    return {
        "selected_pose_index": best_i,
        "score_breakdown": {
            "method": "cluster",
            "n_clusters": cluster_id,
            "largest_cluster_size": cluster_sizes[largest_cluster],
            "largest_cluster_members": members,
            "selected_by_ensemble_within_cluster": True,
        },
        "method_used": "cluster",
    }


def _select_posebusters(
    raw_poses: List[Dict],
    poses_sdf: Path,
    receptor_pdb: Path,
) -> dict:
    """Filter poses by PoseBusters validity, then rank by ensemble score.

    PoseBusters must be installed (``pip install posebusters``). If unavailable,
    falls back to ensemble with a warning.
    """
    try:
        from posebusters import PoseBusters
    except ImportError:
        import warnings
        warnings.warn(
            "select_consensus_pose(method='posebusters'): PoseBusters not installed — falling back to ensemble.",
            stacklevel=3,
        )
        result = _select_ensemble(raw_poses)
        result["method_used"] = "posebusters→ensemble(no_posebusters)"
        return result

    # Run PoseBusters on all poses in the SDF
    try:
        pb = PoseBusters(mode="dock")
        pb_results = pb.bust(str(poses_sdf), str(receptor_pdb))
        # pb_results is a DataFrame with a boolean 'pb_valid' column
        valid_mask = pb_results["pb_valid"].tolist()
    except Exception as exc:
        import warnings
        warnings.warn(
            f"select_consensus_pose(method='posebusters'): PoseBusters run failed ({exc}) — falling back to ensemble.",
            stacklevel=3,
        )
        result = _select_ensemble(raw_poses)
        result["method_used"] = "posebusters→ensemble(pb_error)"
        return result

    # Align mask length with poses list
    valid_indices = [
        i for i, v in enumerate(valid_mask[:len(raw_poses)]) if v
    ]
    if not valid_indices:
        import warnings
        warnings.warn(
            "select_consensus_pose(method='posebusters'): no PB-valid poses found — falling back to ensemble.",
            stacklevel=3,
        )
        result = _select_ensemble(raw_poses)
        result["method_used"] = "posebusters→ensemble(no_valid_pose)"
        return result

    best_i = max(valid_indices, key=lambda i: _ensemble_score(raw_poses[i]))
    return {
        "selected_pose_index": best_i,
        "score_breakdown": {
            "method": "posebusters",
            "n_pb_valid": len(valid_indices),
            "valid_indices": valid_indices,
        },
        "method_used": "posebusters",
    }


def _select_plip(
    raw_poses: List[Dict],
    poses_sdf: Path,
    receptor_pdb: Path,
) -> dict:
    """Rank poses by PLIP interaction count; higher contact count is preferred.

    PLIP must be installed (``pip install plip``). If unavailable,
    falls back to ensemble with a warning.
    """
    try:
        from plip.structure.preparation import PDBComplex
    except ImportError:
        import warnings
        warnings.warn(
            "select_consensus_pose(method='plip'): PLIP not installed — falling back to ensemble.",
            stacklevel=3,
        )
        result = _select_ensemble(raw_poses)
        result["method_used"] = "plip→ensemble(no_plip)"
        return result

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        import warnings
        warnings.warn(
            "select_consensus_pose(method='plip'): RDKit not available — falling back to ensemble.",
            stacklevel=3,
        )
        result = _select_ensemble(raw_poses)
        result["method_used"] = "plip→ensemble(no_rdkit)"
        return result

    import tempfile, os

    receptor_text = Path(receptor_pdb).read_text(encoding="utf-8", errors="replace")
    suppl = Chem.SDMolSupplier(str(poses_sdf), removeHs=False, sanitize=False)
    mols = [m for m in suppl if m is not None]

    contact_counts = []
    for i, mol in enumerate(mols[:len(raw_poses)]):
        try:
            with tempfile.NamedTemporaryFile(
                suffix=".pdb", mode="w", delete=False, encoding="utf-8"
            ) as tmp:
                lig_pdb = Chem.MolToPDBBlock(mol) or ""
                tmp.write(receptor_text)
                tmp.write("\n")
                tmp.write(lig_pdb)
                tmp_path = tmp.name

            cmplx = PDBComplex()
            cmplx.load_pdb(tmp_path)
            cmplx.analyze()
            # count total interactions across all binding sites
            n_contacts = sum(
                len(bs.all_itypes) for bs in cmplx.interaction_sets.values()
            )
        except Exception:
            n_contacts = 0
        finally:
            try:
                os.unlink(tmp_path)
            except Exception:
                pass
        contact_counts.append(n_contacts)

    if not contact_counts or max(contact_counts) == 0:
        import warnings
        warnings.warn(
            "select_consensus_pose(method='plip'): no interactions detected — falling back to ensemble.",
            stacklevel=3,
        )
        result = _select_ensemble(raw_poses)
        result["method_used"] = "plip→ensemble(no_contacts)"
        return result

    best_i = int(max(range(len(contact_counts)), key=lambda i: contact_counts[i]))
    return {
        "selected_pose_index": best_i,
        "score_breakdown": {
            "method": "plip",
            "contact_counts": contact_counts,
            "max_contacts": contact_counts[best_i],
        },
        "method_used": "plip",
    }


def _select_consensus(
    raw_poses: List[Dict],
    poses_sdf: Path,
    receptor_pdb: Path,
    rmsd_cluster_threshold: float = 2.0,
) -> dict:
    """Consensus vote of ensemble + cluster + posebusters.

    If all three methods agree on the same pose, that pose is returned.
    If they disagree, falls back to ensemble.
    """
    r_ensemble   = _select_ensemble(raw_poses)
    r_cluster    = _select_cluster(raw_poses, poses_sdf, rmsd_cluster_threshold)
    r_pb         = _select_posebusters(raw_poses, poses_sdf, receptor_pdb)

    i_ens  = r_ensemble["selected_pose_index"]
    i_clu  = r_cluster["selected_pose_index"]
    i_pb   = r_pb["selected_pose_index"]

    votes = [i_ens, i_clu, i_pb]
    from collections import Counter
    vote_counts = Counter(votes)
    winner, count = vote_counts.most_common(1)[0]

    if count >= 2:
        # majority (≥2 of 3) agree
        return {
            "selected_pose_index": winner,
            "score_breakdown": {
                "method": "consensus",
                "votes": {"ensemble": i_ens, "cluster": i_clu, "posebusters": i_pb},
                "agreement": count,
            },
            "method_used": "consensus",
        }
    else:
        # full disagreement → ensemble fallback
        import warnings
        warnings.warn(
            "select_consensus_pose(method='consensus'): all three methods disagree — falling back to ensemble.",
            stacklevel=3,
        )
        result = _select_ensemble(raw_poses)
        result["score_breakdown"]["votes"] = {
            "ensemble": i_ens, "cluster": i_clu, "posebusters": i_pb,
        }
        result["method_used"] = "consensus→ensemble(disagreement)"
        return result
