"""
dad.core.visualize — Stage 11: Visualization and HTML report generation.

Produces paper-ready outputs per (protein × ligand) pair:
1. Interactive py3Dmol HTML (from GININA_Template.ipynb cells 16-22).
2. ChimeraX .cxc script (from GININA_Template.ipynb cell 23).
3. Project-level HTML report (MultiQC-style summary).

Key visualization parameters (from gnina.ipynb seed notebook)
--------------------------------------------------------------
py3Dmol width   = 900 px
py3Dmol height  = 600 px
cartoon opacity = 0.7
ligand stick radius = 0.25
core residue stick radius = 0.18
contact distance = 5.0 Å (for highlighting residues)

Implementation owner: Mr_Dock (Phase B).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from dad.io import DockingResult, InteractionProfile, RankedPair


# ─────────────────────────────────────────────────────────────────────────────
# py3Dmol visualization
# ─────────────────────────────────────────────────────────────────────────────

def generate_py3dmol_html(
    receptor_pdb: str,
    docked_sdf: str,
    contact_residue_ids: List[int],
    ligand_name: str,
    output_html: str,
    width: int = 900,
    height: int = 600,
    pose_index: int = 0,
) -> str:
    """Generate a standalone HTML file with an interactive py3Dmol view.

    Ports GININA_Template.ipynb cells 20-22:
    - Receptor as spectrum cartoon (opacity 0.7).
    - Best docked pose as green sticks.
    - Contact residues highlighted as gray sticks with labels.
    - Gold cylinders connecting closest atom pairs between ligand and contacts.
    - Ligand label at centroid.

    The output is a self-contained HTML with embedded 3Dmol.js CDN.
    """
    if not Path(receptor_pdb).exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")
    if not Path(docked_sdf).exists():
        raise FileNotFoundError(f"Docked SDF not found: {docked_sdf}")

    receptor_data = Path(receptor_pdb).read_text()
    # extract first model from multi-model SDF
    sdf_content = Path(docked_sdf).read_text()
    models = sdf_content.split("$$$$")
    if pose_index < len(models):
        sdf_first = models[pose_index].strip() + "\n$$$$\n"
    else:
        sdf_first = models[0].strip() + "\n$$$$\n"

    # compute interaction pairs for cylinders (mirrors cell 19-20)
    interaction_pairs = _compute_interaction_pairs(receptor_pdb, docked_sdf, pose_index)
    # compute ligand centroid for label
    from dad.core.interaction import get_sdf_atom_coords
    import numpy as np
    lig_coords = get_sdf_atom_coords(docked_sdf, pose_index)
    centroid = lig_coords.mean(axis=0).tolist() if len(lig_coords) > 0 else [0, 0, 0]

    resi_list = [str(r) for r in contact_residue_ids]
    cylinders_js = _build_cylinders_js(interaction_pairs)

    html = _PY3DMOL_TEMPLATE.format(
        width=width,
        height=height,
        receptor_data=json.dumps(receptor_data),
        sdf_data=json.dumps(sdf_first),
        resi_list=json.dumps(resi_list),
        ligand_name=json.dumps(ligand_name),
        centroid_x=centroid[0],
        centroid_y=centroid[1],
        centroid_z=centroid[2],
        cylinders_js=cylinders_js,
    )

    Path(output_html).parent.mkdir(parents=True, exist_ok=True)
    Path(output_html).write_text(html, encoding="utf-8")
    return output_html


def _compute_interaction_pairs(
    receptor_pdb: str,
    docked_sdf: str,
    pose_index: int = 0,
) -> List[Tuple[List[float], List[float]]]:
    """Compute atom-pair contacts for cylinder drawing (mirrors cell 19)."""
    from dad.core.interaction import get_sdf_atom_coords
    import numpy as np

    lig_coords = get_sdf_atom_coords(docked_sdf, pose_index)
    if len(lig_coords) == 0:
        return []

    try:
        from Bio.PDB import PDBParser
    except ModuleNotFoundError:
        return _compute_interaction_pairs_plain(receptor_pdb, lig_coords)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("r", receptor_pdb)

    pairs = []
    search_distance = 5.0
    for model in structure:
        for chain in model:
            for residue in chain:
                min_dist = float("inf")
                best_l = None
                best_p = None
                for atom in residue:
                    p = atom.get_coord()
                    dists = np.linalg.norm(lig_coords - p, axis=1)
                    idx = int(dists.argmin())
                    d = float(dists[idx])
                    if d < min_dist:
                        min_dist = d
                        best_l = lig_coords[idx].tolist()
                        best_p = p.tolist()
                if min_dist <= search_distance and best_l and best_p:
                    pairs.append((best_l, best_p))
    return pairs


def _compute_interaction_pairs_plain(
    receptor_pdb: str,
    lig_coords,
) -> List[Tuple[List[float], List[float]]]:
    """Fallback atom-pair contacts used when BioPython is unavailable."""
    import numpy as np

    residues: Dict[Tuple[str, int, str], List[List[float]]] = {}
    for line in Path(receptor_pdb).read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            key = (line[21].strip() or "A", int(line[22:26]), line[17:20].strip() or "UNK")
            residues.setdefault(key, []).append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ])
        except ValueError:
            continue

    pairs: List[Tuple[List[float], List[float]]] = []
    for atoms in residues.values():
        min_dist = float("inf")
        best_l = None
        best_p = None
        for p_coord in atoms:
            p = np.array(p_coord, dtype=float)
            dists = np.linalg.norm(lig_coords - p, axis=1)
            idx = int(dists.argmin())
            d = float(dists[idx])
            if d < min_dist:
                min_dist = d
                best_l = lig_coords[idx].tolist()
                best_p = p_coord
        if min_dist <= 5.0 and best_l and best_p:
            pairs.append((best_l, best_p))
    return pairs


def _build_cylinders_js(pairs: List[Tuple[List[float], List[float]]]) -> str:
    lines = []
    for l_coord, p_coord in pairs:
        lines.append(
            f"viewer.addCylinder({{"
            f"start:{{x:{l_coord[0]:.3f},y:{l_coord[1]:.3f},z:{l_coord[2]:.3f}}},"
            f"end:{{x:{p_coord[0]:.3f},y:{p_coord[1]:.3f},z:{p_coord[2]:.3f}}},"
            f"radius:0.03,color:'#FFD700',fromCap:1,toCap:1}});"
        )
    return "\n".join(lines)


_PY3DMOL_TEMPLATE = """\
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>DAD Docking: {ligand_name}</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; background: #fff; }}
  #viewer {{ width: {width}px; height: {height}px; position: relative; }}
  h2 {{ text-align: center; color: #2e7d32; }}
</style>
</head>
<body>
<h2>DAD Docking View: {ligand_name}</h2>
<div id="viewer"></div>
<script>
(function() {{
  var viewer = $3Dmol.createViewer('viewer', {{backgroundColor: 'white'}});

  // receptor
  viewer.addModel({receptor_data}, 'pdb');
  viewer.setStyle({{model: 0}}, {{cartoon: {{color: 'spectrum', opacity: 0.7}}}});

  // docked pose (best / first model)
  viewer.addModel({sdf_data}, 'sdf');
  viewer.setStyle({{model: 1}}, {{stick: {{colorscheme: 'greenCarbon', radius: 0.25}}}});

  // contact residues
  var resiList = {resi_list};
  viewer.addStyle({{model: 0, resi: resiList}},
    {{stick: {{colorscheme: 'lightgrayCarbon', radius: 0.18}}}});
  viewer.addResLabels({{model: 0, resi: resiList}}, {{
    fontSize: 13, fontColor: 'black',
    showBackground: true, backgroundColor: 'white', backgroundOpacity: 0.8
  }});

  // ligand label at centroid
  viewer.addLabel({ligand_name}, {{
    position: {{x: {centroid_x}, y: {centroid_y}, z: {centroid_z}}},
    backgroundColor: '#2e7d32', fontColor: 'white', fontSize: 14,
    backgroundOpacity: 0.9, alignment: 'center'
  }});

  // interaction cylinders
  {cylinders_js}

  viewer.zoomTo({{model: 1}});
  viewer.render();
}})();
</script>
</body>
</html>
"""


# ─────────────────────────────────────────────────────────────────────────────
# ChimeraX .cxc script generation
# ─────────────────────────────────────────────────────────────────────────────

def generate_chimerax_cxc(
    receptor_pdb: str,
    docked_sdf: str,
    contact_residue_ids: List[int],
    ligand_name: str,
    output_cxc: str,
    complex_pdb: Optional[str] = None,
) -> str:
    """Generate a ChimeraX .cxc script for paper-quality 3D visualization.

    Ports GININA_Template.ipynb cell 23 (save_complex_with_style):
    - Merge receptor + ligand if complex_pdb is not provided.
    - Write .cxc mirroring the AW1_ref PDB_file/*.cxc style:
        open <complex.pdb>
        hide all
        show /A cartoon
        color bychain; transparency 30 target c
        show ligand sticks; color ligand green
        show /<res_list> sticks; color /<res_list> lightgray
        label commands for residues and ligand
        set bgcolor white
        view ligand
    """
    if not Path(receptor_pdb).exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")
    if not Path(docked_sdf).exists():
        raise FileNotFoundError(f"Docked SDF not found: {docked_sdf}")

    if complex_pdb is None:
        output_dir = Path(output_cxc).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        complex_pdb = str(output_dir / (Path(output_cxc).stem + "_complex.pdb"))
        merge_receptor_ligand_pdb(receptor_pdb, docked_sdf, complex_pdb)

    res_list_str = ",".join(str(r) for r in contact_residue_ids)
    complex_name = Path(complex_pdb).name

    lines = [
        f"open {complex_name}",
        "hide all",
        "show /A cartoon",
        "color bychain",
        "transparency 30 target c",
        "show ligand sticks",
        "color ligand green",
        f"label ligand text '{ligand_name}' size 14 color green",
    ]

    if res_list_str:
        lines += [
            f"show /{res_list_str} sticks",
            f"color /{res_list_str} lightgray",
            f"label /{res_list_str} text '{{0.name}}{{0.number}}' size 13 color black",
        ]

    lines += [
        "set bgcolor white",
        "view ligand",
    ]

    Path(output_cxc).parent.mkdir(parents=True, exist_ok=True)
    Path(output_cxc).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_cxc


def merge_receptor_ligand_pdb(
    receptor_pdb: str,
    docked_sdf: str,
    output_pdb: str,
    pose_index: int = 0,
) -> str:
    """Merge receptor PDB and docked ligand SDF into a single PDB for ChimeraX.

    Steps from GININA_Template.ipynb cell 23:
    1. Read receptor PDB, strip END/CONECT records.
    2. Read ligand SDF (best pose), convert to PDB block via Chem.MolToPDBBlock.
    3. Concatenate and write to output_pdb.
    """
    if not Path(receptor_pdb).exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")
    if not Path(docked_sdf).exists():
        raise FileNotFoundError(f"Docked SDF not found: {docked_sdf}")

    try:
        from rdkit import Chem
    except ModuleNotFoundError:
        return _merge_receptor_ligand_pdb_plain(
            receptor_pdb=receptor_pdb,
            docked_sdf=docked_sdf,
            output_pdb=output_pdb,
            pose_index=pose_index,
        )

    with open(receptor_pdb, "r") as fh:
        rec_lines = [l for l in fh.readlines() if not l.rstrip().startswith(("END", "CONECT"))]

    sdf_content = Path(docked_sdf).read_text()
    models = sdf_content.split("$$$$")
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

    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as fh:
        fh.writelines(rec_lines)
        fh.writelines(lig_lines)
        fh.write("END\n")

    return output_pdb


def _merge_receptor_ligand_pdb_plain(
    receptor_pdb: str,
    docked_sdf: str,
    output_pdb: str,
    pose_index: int = 0,
) -> str:
    receptor_lines = [
        line for line in Path(receptor_pdb).read_text(errors="replace").splitlines()
        if not line.startswith(("END", "CONECT"))
    ]
    max_serial = 0
    for line in receptor_lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                max_serial = max(max_serial, int(line[6:11]))
            except ValueError:
                pass

    ligand_atoms, ligand_bonds = _parse_sdf_atoms_bonds(docked_sdf, pose_index)
    lines = list(receptor_lines)
    serial_map: Dict[int, int] = {}
    for idx, atom in enumerate(ligand_atoms, start=1):
        serial = max_serial + idx
        serial_map[idx] = serial
        elem = (atom["element"] or "C").upper()[:2]
        atom_name = f"{elem}{idx}"[:4]
        lines.append(
            f"HETATM{serial:5d} {atom_name:<4s} LIG Z   1    "
            f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}"
            f"  1.00 20.00          {elem:>2s}"
        )

    for atom_a, atom_b in ligand_bonds:
        if atom_a in serial_map and atom_b in serial_map:
            lines.append(f"CONECT{serial_map[atom_a]:5d}{serial_map[atom_b]:5d}")

    lines.append("END")
    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    Path(output_pdb).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_pdb


def _parse_sdf_atoms_bonds(
    sdf_path: str,
    pose_index: int = 0,
) -> Tuple[List[Dict[str, object]], List[Tuple[int, int]]]:
    text = Path(sdf_path).read_text(errors="replace")
    models = [m for m in text.split("$$$$") if m.strip()]
    if not models:
        return [], []
    pose_index = max(0, min(pose_index, len(models) - 1))
    lines = models[pose_index].splitlines()

    counts_idx = None
    atom_count = 0
    bond_count = 0
    for i, line in enumerate(lines[:20]):
        if "V2000" not in line:
            continue
        parts = line.split()
        try:
            atom_count = int(parts[0])
            bond_count = int(parts[1])
            counts_idx = i
        except (ValueError, IndexError):
            pass
        break
    if counts_idx is None or atom_count <= 0:
        return [], []

    atoms: List[Dict[str, object]] = []
    atom_start = counts_idx + 1
    for line in lines[atom_start:atom_start + atom_count]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            atoms.append({
                "x": float(parts[0]),
                "y": float(parts[1]),
                "z": float(parts[2]),
                "element": "".join(ch for ch in parts[3] if ch.isalpha())[:2] or "C",
            })
        except ValueError:
            continue

    bonds: List[Tuple[int, int]] = []
    bond_start = atom_start + atom_count
    for line in lines[bond_start:bond_start + bond_count]:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            bonds.append((int(parts[0]), int(parts[1])))
        except ValueError:
            continue
    return atoms, bonds


# ─────────────────────────────────────────────────────────────────────────────
# Project-level HTML report
# ─────────────────────────────────────────────────────────────────────────────

def generate_html_report(
    ranked_pairs: List[RankedPair],
    interaction_profiles: List[InteractionProfile],
    html_dir: str,
    output_html: str = "results/report/dad_report.html",
    top_k: int = 10,
) -> str:
    """Generate a project-level HTML summary report (MultiQC-style).

    Contents:
    - Run metadata (date, n_proteins, n_ligands).
    - Top-K ranked pairs table with links to per-pair HTML.
    - Per-pair py3Dmol viewer iframes for top-K.
    """
    import datetime

    profile_map = {
        (p.seq_id, p.lig_id, p.pocket_rank): p
        for p in interaction_profiles
    }

    proteins = sorted({rp.seq_id for rp in ranked_pairs})
    ligands = sorted({rp.lig_id for rp in ranked_pairs})
    top = ranked_pairs[:top_k]

    rows_html = []
    for rp in top:
        prof = profile_map.get((rp.seq_id, rp.lig_id, rp.pocket_rank))
        pair_html = Path(html_dir) / f"{rp.seq_id}__{rp.lig_id}__pocket{rp.pocket_rank}.html"
        link = pair_html.name if pair_html.exists() else "#"
        rows_html.append(
            f"<tr>"
            f"<td>{rp.overall_rank}</td>"
            f"<td><a href='{link}'>{rp.seq_id}</a></td>"
            f"<td>{rp.lig_id}</td>"
            f"<td>{rp.pocket_rank}</td>"
            f"<td>{rp.vina_affinity:.3f}</td>"
            f"<td>{rp.cnn_pose_score:.4f}</td>"
            f"<td>{rp.cnn_affinity:.3f}</td>"
            f"<td>{rp.composite_rank_score:.4f}</td>"
            f"<td>{rp.n_contacts if rp.n_contacts is not None else '-'}</td>"
            f"<td>{rp.n_hbonds if rp.n_hbonds is not None else '-'}</td>"
            f"</tr>"
        )

    iframes_html = []
    for rp in top:
        pair_html = Path(html_dir) / f"{rp.seq_id}__{rp.lig_id}__pocket{rp.pocket_rank}.html"
        if pair_html.exists():
            iframes_html.append(
                f"<h3>{rp.seq_id} × {rp.lig_id} (pocket {rp.pocket_rank})</h3>"
                f"<iframe src='{pair_html.name}' width='920' height='640' frameborder='0'></iframe>"
            )

    report_html = _REPORT_TEMPLATE.format(
        date=datetime.date.today().isoformat(),
        n_proteins=len(proteins),
        n_ligands=len(ligands),
        n_pairs=len(ranked_pairs),
        table_rows="\n".join(rows_html),
        iframes="\n".join(iframes_html),
    )

    Path(output_html).parent.mkdir(parents=True, exist_ok=True)
    Path(output_html).write_text(report_html, encoding="utf-8")
    return output_html


_REPORT_TEMPLATE = """\
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>DAD Docking Report</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 20px; background: #f9f9f9; }}
  h1 {{ color: #1a237e; }}
  h2 {{ color: #2e7d32; border-bottom: 1px solid #ccc; }}
  table {{ border-collapse: collapse; width: 100%; background: #fff; }}
  th, td {{ border: 1px solid #ddd; padding: 6px 10px; text-align: right; }}
  th {{ background: #1a237e; color: #fff; text-align: center; }}
  td:nth-child(2), td:nth-child(3) {{ text-align: left; }}
  tr:nth-child(even) {{ background: #f5f5f5; }}
  .meta {{ font-size: 0.9em; color: #555; margin-bottom: 20px; }}
  iframe {{ border: 1px solid #ccc; margin-bottom: 20px; }}
</style>
</head>
<body>
<h1>DAD Dynamic Affinity Dock — Report</h1>
<div class="meta">
  Generated: {date} &nbsp;|&nbsp;
  Proteins: {n_proteins} &nbsp;|&nbsp;
  Ligands: {n_ligands} &nbsp;|&nbsp;
  Pairs evaluated: {n_pairs}
</div>
<h2>Top-Ranked Pairs</h2>
<table>
<tr>
  <th>Rank</th><th>Protein</th><th>Ligand</th><th>Pocket</th>
  <th>Vina (kcal/mol)</th><th>CNN Pose</th><th>CNN Affinity</th>
  <th>Composite Score</th><th>Contacts</th><th>H-bonds</th>
</tr>
{table_rows}
</table>
<h2>3D Views (Top Pairs)</h2>
{iframes}
</body>
</html>
"""


# ─────────────────────────────────────────────────────────────────────────────
# Batch visualization entry point (called by Snakemake rule)
# ─────────────────────────────────────────────────────────────────────────────

def run_visualization_batch(
    docking_results: List[DockingResult],
    interaction_profiles: List[InteractionProfile],
    ranked_pairs: List[RankedPair],
    structures: Dict,
    config: Dict,
    output_dir: str = "results/visualization",
) -> Dict[str, List[str]]:
    """Run visualization for all (protein × ligand) pairs.

    Orchestrates generate_py3dmol_html + generate_chimerax_cxc for every
    DockingResult, then calls generate_html_report.
    """
    vis_cfg = config.get("visualization", {})
    width = vis_cfg.get("width", 900)
    height = vis_cfg.get("height", 600)

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    profile_map = {
        (p.seq_id, p.lig_id, p.pocket_rank): p
        for p in interaction_profiles
    }

    html_files: List[str] = []
    cxc_files: List[str] = []

    for result in docking_results:
        struct = structures.get(result.seq_id)
        if struct is None:
            continue
        receptor_pdb = getattr(struct, "pdb_path", str(struct))
        profile = profile_map.get((result.seq_id, result.lig_id, result.pocket_rank))
        contact_ids = (
            [c.res_id for c in profile.contact_residues] if profile else []
        )

        stem = f"{result.seq_id}__{result.lig_id}__pocket{result.pocket_rank}"
        html_path = str(out_path / f"{stem}.html")
        cxc_path = str(out_path / f"{stem}.cxc")

        try:
            generate_py3dmol_html(
                receptor_pdb=receptor_pdb,
                docked_sdf=result.output_sdf,
                contact_residue_ids=contact_ids,
                ligand_name=result.lig_id,
                output_html=html_path,
                width=width,
                height=height,
            )
            html_files.append(html_path)
        except Exception as exc:
            import logging
            logging.getLogger(__name__).warning("py3dmol failed for %s: %s", stem, exc)

        try:
            generate_chimerax_cxc(
                receptor_pdb=receptor_pdb,
                docked_sdf=result.output_sdf,
                contact_residue_ids=contact_ids,
                ligand_name=result.lig_id,
                output_cxc=cxc_path,
            )
            cxc_files.append(cxc_path)
        except Exception as exc:
            import logging
            logging.getLogger(__name__).warning("chimerax cxc failed for %s: %s", stem, exc)

    report_html = str(out_path / "dad_report.html")
    try:
        generate_html_report(
            ranked_pairs=ranked_pairs,
            interaction_profiles=interaction_profiles,
            html_dir=str(out_path),
            output_html=report_html,
        )
    except Exception as exc:
        import logging
        logging.getLogger(__name__).warning("html report generation failed: %s", exc)
        report_html = ""

    return {
        "html_files": html_files,
        "cxc_files": cxc_files,
        "report_html": report_html,
    }
