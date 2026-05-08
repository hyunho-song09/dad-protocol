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
# Multi-pose animation HTML (GININA_Template cells 16-17)
# ─────────────────────────────────────────────────────────────────────────────

def generate_multi_pose_animation_html(
    receptor_pdb: str,
    docked_sdf: str,
    output_html: str,
    width: int = 900,
    height: int = 600,
    interval_ms: int = 1200,
    cartoon_style: str = "spectrum",
) -> str:
    """Generate py3Dmol HTML with all GNINA poses as animation frames.

    Ports GININA_Template.ipynb cells 16-17 (addModelsAsFrames + animate).
    The receptor is shown as a cartoon and each docked pose cycles as a frame.
    Complements generate_py3dmol_html (single-pose static) — both can be used
    together per case.
    """
    if not Path(receptor_pdb).exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")
    if not Path(docked_sdf).exists():
        raise FileNotFoundError(f"Docked SDF not found: {docked_sdf}")

    receptor_data = Path(receptor_pdb).read_text(encoding="utf-8", errors="replace")
    sdf_data = Path(docked_sdf).read_text(encoding="utf-8", errors="replace")

    # count poses for display
    n_poses = len([m for m in sdf_data.split("$$$$") if m.strip()])

    cartoon_color = json.dumps(cartoon_style)

    html = _ANIMATION_TEMPLATE.format(
        width=width,
        height=height,
        receptor_data=json.dumps(receptor_data),
        sdf_data=json.dumps(sdf_data),
        interval_ms=interval_ms,
        cartoon_color=cartoon_color,
        n_poses=n_poses,
        docked_sdf_name=Path(docked_sdf).name,
    )

    Path(output_html).parent.mkdir(parents=True, exist_ok=True)
    Path(output_html).write_text(html, encoding="utf-8")
    return output_html


_ANIMATION_TEMPLATE = """\
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>DAD Multi-Pose Animation ({n_poses} poses)</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; background: #fff; margin: 12px; }}
  #viewer {{ width: {width}px; height: {height}px; position: relative; }}
  h3 {{ color: #1a237e; margin-bottom: 4px; }}
  .info {{ font-size: 0.85em; color: #555; margin-bottom: 8px; }}
  button {{ margin: 4px 2px; padding: 4px 12px; cursor: pointer; }}
</style>
</head>
<body>
<h3>DAD Multi-Pose Animation</h3>
<div class="info">
  {n_poses} GNINA poses from <code>{docked_sdf_name}</code> — cycling every {interval_ms} ms
</div>
<div id="viewer"></div>
<div>
  <button onclick="animViewer.animate({{interval:{interval_ms}}})">Play</button>
  <button onclick="animViewer.stopAnimate()">Pause</button>
  <button onclick="animViewer.setFrame(0); animViewer.render();">Reset</button>
</div>
<script>
(function() {{
  var viewer = $3Dmol.createViewer('viewer', {{backgroundColor: 'white'}});
  window.animViewer = viewer;

  // receptor (model 0)
  viewer.addModel({receptor_data}, 'pdb');
  viewer.setStyle({{model: 0}}, {{
    cartoon: {{color: {cartoon_color}, opacity: 0.70}},
    stick: {{radius: 0.10, colorscheme: 'grayCarbon'}}
  }});

  // all docked poses as animation frames (model 1+)
  viewer.addModelsAsFrames({sdf_data}, 'sdf');
  viewer.setStyle({{model: -1}}, {{stick: {{colorscheme: 'greenCarbon', radius: 0.22}}}});

  viewer.zoomTo({{model: -1}});
  viewer.rotate(90, 'y');
  viewer.animate({{interval: {interval_ms}, loop: 'forward'}});
  viewer.render();
}})();
</script>
</body>
</html>
"""


# ─────────────────────────────────────────────────────────────────────────────
# Multi-case grid HTML
# ─────────────────────────────────────────────────────────────────────────────

def generate_multi_case_grid_html(
    cases: List[Dict],
    output_html: str,
    columns: int = 3,
    width: int = 350,
    height: int = 280,
) -> str:
    """Composite HTML showing all selected (protein × ligand) pairs in a grid.

    Each cell contains a small py3Dmol view with the best docked pose and a
    score summary text block below it. Self-contained with 3Dmol.js CDN.

    Each entry in ``cases`` must have keys:
        case_id, receptor_pdb, docked_sdf, ligand_name
    Optional keys: contact_residues (List[int]), score_summary (str).
    """
    cell_blocks = []
    for case in cases:
        case_id = case.get("case_id", "case")
        receptor_pdb = case.get("receptor_pdb", "")
        docked_sdf = case.get("docked_sdf", "")
        ligand_name = case.get("ligand_name", case_id)
        score_summary = case.get("score_summary", "")
        contact_ids = case.get("contact_residues", [])

        # Only embed data for cases where files exist
        if Path(receptor_pdb).exists() and Path(docked_sdf).exists():
            rec_data = json.dumps(
                Path(receptor_pdb).read_text(encoding="utf-8", errors="replace")
            )
            sdf_content = Path(docked_sdf).read_text(encoding="utf-8", errors="replace")
            models = sdf_content.split("$$$$")
            first_sdf = (models[0].strip() + "\n$$$$\n") if models else ""
            sdf_data = json.dumps(first_sdf)
            resi_list = json.dumps([str(r) for r in contact_ids])
            viewer_id = f"v_{case_id.replace('-', '_').replace('.', '_')}"

            viewer_js = f"""
<div id="{viewer_id}" style="width:{width}px;height:{height}px;position:relative;"></div>
<script>
(function() {{
  var viewer = $3Dmol.createViewer('{viewer_id}', {{backgroundColor:'white'}});
  viewer.addModel({rec_data}, 'pdb');
  viewer.setStyle({{model:0}}, {{cartoon:{{color:'spectrum',opacity:0.75}}}});
  viewer.addModel({sdf_data}, 'sdf');
  viewer.setStyle({{model:1}}, {{stick:{{colorscheme:'greenCarbon',radius:0.22}}}});
  var resi = {resi_list};
  if (resi.length > 0) {{
    viewer.addStyle({{model:0,resi:resi}}, {{stick:{{colorscheme:'lightgrayCarbon',radius:0.16}}}});
  }}
  viewer.zoomTo({{model:1}});
  viewer.render();
}})();
</script>"""
        else:
            viewer_js = (
                f'<div style="width:{width}px;height:{height}px;background:#f5f5f5;'
                f'display:flex;align-items:center;justify-content:center;'
                f'font-size:12px;color:#888;">'
                f'Files not found</div>'
            )

        score_html = (
            f'<div style="font-size:11px;color:#333;padding:4px;'
            f'max-width:{width}px;overflow:hidden;white-space:pre-wrap;">'
            f'{score_summary}</div>'
        ) if score_summary else ""

        cell_blocks.append(
            f'<div style="display:inline-block;vertical-align:top;'
            f'margin:8px;border:1px solid #ddd;padding:6px;background:#fff;">'
            f'<div style="font-weight:bold;font-size:12px;color:#1a237e;'
            f'margin-bottom:4px;">{case_id}</div>'
            f'{viewer_js}'
            f'{score_html}'
            f'</div>'
        )

    grid_html = _GRID_TEMPLATE.format(
        n_cases=len(cases),
        columns=columns,
        cells="\n".join(cell_blocks),
    )

    Path(output_html).parent.mkdir(parents=True, exist_ok=True)
    Path(output_html).write_text(grid_html, encoding="utf-8")
    return output_html


_GRID_TEMPLATE = """\
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>DAD Docking Grid ({n_cases} cases)</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; background: #f9f9f9; margin: 16px; }}
  h2 {{ color: #1a237e; }}
  .grid {{ display: flex; flex-wrap: wrap; gap: 0; }}
</style>
</head>
<body>
<h2>DAD Docking Overview — {n_cases} pairs</h2>
<div class="grid">
{cells}
</div>
</body>
</html>
"""


# ─────────────────────────────────────────────────────────────────────────────
# Static PNG export
# ─────────────────────────────────────────────────────────────────────────────

def export_pose_static_png(
    receptor_pdb: str,
    docked_sdf: str,
    output_png: str,
    contact_residue_ids: List[int],
    ligand_name: str,
    width: int = 900,
    height: int = 600,
    pose_index: int = 0,
) -> str:
    """Generate a static PNG of the docked pose for paper figures.

    Strategy (in order):
    1. py3Dmol ``view.png()`` — works in Colab with Chrome headless.
    2. matplotlib 2D contact-distance heatmap fallback (no 3D render needed).
    3. Last resort: write a companion HTML and return its path with a note.
    """
    if not Path(receptor_pdb).exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")
    if not Path(docked_sdf).exists():
        raise FileNotFoundError(f"Docked SDF not found: {docked_sdf}")

    Path(output_png).parent.mkdir(parents=True, exist_ok=True)

    # ── Strategy 1: py3Dmol headless render ───────────────────────────────────
    try:
        import py3Dmol as _py3dmol
        from dad.core.interaction import get_sdf_atom_coords
        import numpy as np

        rec_data = Path(receptor_pdb).read_text(encoding="utf-8", errors="replace")
        sdf_content = Path(docked_sdf).read_text(encoding="utf-8", errors="replace")
        models = sdf_content.split("$$$$")
        idx = min(pose_index, len(models) - 1)
        sdf_first = models[idx].strip() + "\n$$$$\n"

        view = _py3dmol.view(width=width, height=height)
        view.addModel(rec_data, "pdb")
        view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum", "opacity": 0.7}})
        view.addModel(sdf_first, "sdf")
        view.setStyle({"model": 1}, {"stick": {"colorscheme": "greenCarbon", "radius": 0.25}})
        resi = [str(r) for r in contact_residue_ids]
        if resi:
            view.addStyle({"model": 0, "resi": resi},
                          {"stick": {"colorscheme": "lightgrayCarbon", "radius": 0.18}})
        view.zoomTo({"model": 1})
        view.rotate(90, "y")

        if hasattr(view, "png"):
            png_data = view.png()
            if png_data:
                import base64
                raw = png_data.split(",", 1)[-1]
                Path(output_png).write_bytes(base64.b64decode(raw))
                return output_png
    except Exception:
        pass

    # ── Strategy 2: matplotlib contact-distance heatmap fallback ─────────────
    try:
        return _export_contact_heatmap_png(
            receptor_pdb=receptor_pdb,
            docked_sdf=docked_sdf,
            output_png=output_png,
            contact_residue_ids=contact_residue_ids,
            ligand_name=ligand_name,
            pose_index=pose_index,
            width=width,
            height=height,
        )
    except Exception:
        pass

    # ── Strategy 3: companion HTML fallback ───────────────────────────────────
    import warnings
    companion_html = str(Path(output_png).with_suffix(".html"))
    warnings.warn(
        f"export_pose_static_png: neither py3Dmol-png nor matplotlib succeeded. "
        f"Writing companion HTML to {companion_html}. "
        f"Open in Chrome/Colab for interactive view.",
        stacklevel=2,
    )
    generate_py3dmol_html(
        receptor_pdb=receptor_pdb,
        docked_sdf=docked_sdf,
        contact_residue_ids=contact_residue_ids,
        ligand_name=ligand_name,
        output_html=companion_html,
        width=width,
        height=height,
        pose_index=pose_index,
    )
    # write a placeholder PNG with note
    _write_placeholder_png(output_png, ligand_name, companion_html)
    return output_png


def _export_contact_heatmap_png(
    receptor_pdb: str,
    docked_sdf: str,
    output_png: str,
    contact_residue_ids: List[int],
    ligand_name: str,
    pose_index: int = 0,
    width: int = 900,
    height: int = 600,
) -> str:
    """Write a matplotlib contact-distance heatmap as PNG fallback."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from dad.core.interaction import get_sdf_atom_coords

    lig_coords = get_sdf_atom_coords(docked_sdf, pose_index)
    if len(lig_coords) == 0:
        raise ValueError("No ligand coordinates")

    # collect contact residue atoms from PDB
    contact_set = set(str(r) for r in contact_residue_ids)
    residue_atoms: Dict[str, List] = {}
    for line in Path(receptor_pdb).read_text(errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            res_id = line[22:26].strip()
            if res_id not in contact_set:
                continue
            res_name = line[17:20].strip()
            label = f"{res_name}{res_id}"
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            residue_atoms.setdefault(label, []).append([x, y, z])
        except (ValueError, IndexError):
            continue

    if not residue_atoms:
        # nothing to plot — create minimal figure
        labels = [str(r) for r in contact_residue_ids[:20]]
        residue_atoms = {lbl: [[0, 0, 0]] for lbl in labels}

    res_labels = sorted(residue_atoms.keys())
    n_res = len(res_labels)
    n_lig = len(lig_coords)

    # distance matrix: res × ligand_atom
    dist_matrix = np.zeros((n_res, n_lig))
    for i, label in enumerate(res_labels):
        atoms = np.array(residue_atoms[label])
        for j, lc in enumerate(lig_coords):
            dist_matrix[i, j] = float(np.min(np.linalg.norm(atoms - lc, axis=1)))

    dpi = 100
    fig_w = max(6, min(n_lig * 0.4, width / dpi))
    fig_h = max(4, min(n_res * 0.35, height / dpi))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    im = ax.imshow(dist_matrix, aspect="auto", cmap="RdYlGn_r", vmin=0, vmax=8)
    ax.set_xticks(range(n_lig))
    ax.set_xticklabels([f"A{j+1}" for j in range(n_lig)], fontsize=7, rotation=60)
    ax.set_yticks(range(n_res))
    ax.set_yticklabels(res_labels, fontsize=8)
    ax.set_xlabel("Ligand atom", fontsize=9)
    ax.set_ylabel("Contact residue", fontsize=9)
    ax.set_title(f"Contact distances — {ligand_name} (pose {pose_index})", fontsize=10)
    fig.colorbar(im, ax=ax, label="Distance (Å)")
    fig.tight_layout()
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output_png


def _write_placeholder_png(output_png: str, ligand_name: str, companion_html: str) -> None:
    """Write a minimal white-with-text PNG when no render is available."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4, 2), dpi=72)
        ax.axis("off")
        ax.text(0.5, 0.6, f"PNG export pending — {ligand_name}",
                ha="center", va="center", fontsize=10)
        ax.text(0.5, 0.35, f"Open {Path(companion_html).name} in browser",
                ha="center", va="center", fontsize=8, color="#555")
        fig.savefig(output_png, dpi=72, bbox_inches="tight")
        plt.close(fig)
    except Exception:
        Path(output_png).write_bytes(b"")


# ─────────────────────────────────────────────────────────────────────────────
# 2D interaction map (RDKit + matplotlib)
# ─────────────────────────────────────────────────────────────────────────────

def generate_2d_interaction_map(
    docked_sdf: str,
    contact_residue_labels: List[str],
    contact_distances: List[float],
    output_png: str,
    pose_index: int = 0,
    width: int = 800,
    height: int = 600,
) -> str:
    """RDKit 2D ligand depiction + contact residue bar chart.

    Left panel: RDKit Chem.Draw 2D ligand structure.
    Right panel: horizontal bar chart of residue labels sorted by distance.
    Bars coloured by distance (green < 3 Å, yellow 3-4 Å, orange > 4 Å).

    Falls back to a text-only PNG when RDKit or matplotlib is unavailable.
    """
    if not Path(docked_sdf).exists():
        raise FileNotFoundError(f"Docked SDF not found: {docked_sdf}")

    Path(output_png).parent.mkdir(parents=True, exist_ok=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        import warnings
        warnings.warn(
            "generate_2d_interaction_map: matplotlib not available — writing text fallback.",
            stacklevel=2,
        )
        txt_path = str(Path(output_png).with_suffix(".txt"))
        Path(txt_path).write_text(
            "\n".join(f"{lbl}: {d:.2f} Å" for lbl, d in
                      zip(contact_residue_labels, contact_distances)),
            encoding="utf-8",
        )
        return txt_path

    # ── RDKit 2D ligand depiction ─────────────────────────────────────────────
    lig_img = None
    try:
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        from PIL import Image
        import io

        sdf_text = Path(docked_sdf).read_text(encoding="utf-8", errors="replace")
        models = sdf_text.split("$$$$")
        idx = min(pose_index, len(models) - 1)
        mol_block = models[idx].strip() + "\n$$$$\n"
        suppl = Chem.SDMolSupplier()
        suppl.SetData(mol_block)
        mol = next((m for m in suppl if m is not None), None)
        if mol is not None:
            mol = Chem.RemoveHs(mol)
            draw_w, draw_h = 380, height - 40
            drawer = rdMolDraw2D.MolDraw2DSVG(draw_w, draw_h)
            drawer.drawOptions().addStereoAnnotation = True
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            # convert SVG → PIL image via cairosvg if available, else skip
            try:
                import cairosvg
                png_bytes = cairosvg.svg2png(bytestring=svg.encode(), output_width=draw_w)
                lig_img = Image.open(io.BytesIO(png_bytes))
            except Exception:
                lig_img = None
    except Exception:
        lig_img = None

    # ── Build figure ──────────────────────────────────────────────────────────
    dpi = 100
    fig_w = width / dpi
    fig_h = height / dpi

    if lig_img is not None:
        fig, (ax_lig, ax_bar) = plt.subplots(
            1, 2, figsize=(fig_w, fig_h), dpi=dpi,
            gridspec_kw={"width_ratios": [2, 3]}
        )
        ax_lig.imshow(lig_img)
        ax_lig.axis("off")
        ax_lig.set_title("2D Structure", fontsize=9)
    else:
        fig, ax_bar = plt.subplots(1, 1, figsize=(fig_w * 0.65, fig_h), dpi=dpi)

    # ── Contact bar chart ─────────────────────────────────────────────────────
    if contact_residue_labels and contact_distances:
        paired = sorted(
            zip(contact_residue_labels, contact_distances),
            key=lambda x: x[1],
        )
        labels_sorted = [p[0] for p in paired]
        dists_sorted = [p[1] for p in paired]
        colors = [
            "#2e7d32" if d < 3.0 else ("#f9a825" if d < 4.0 else "#e64a19")
            for d in dists_sorted
        ]
        y_pos = list(range(len(labels_sorted)))
        ax_bar.barh(y_pos, dists_sorted, color=colors, height=0.6, edgecolor="white")
        ax_bar.set_yticks(y_pos)
        ax_bar.set_yticklabels(labels_sorted, fontsize=8)
        ax_bar.set_xlabel("Min. distance to ligand (Å)", fontsize=9)
        ax_bar.axvline(x=3.0, color="#2e7d32", linestyle="--", linewidth=0.8, alpha=0.7)
        ax_bar.axvline(x=4.0, color="#f9a825", linestyle="--", linewidth=0.8, alpha=0.7)
        ax_bar.set_xlim(0, max(dists_sorted) * 1.1 + 0.5)
        ax_bar.invert_yaxis()
        ax_bar.set_title("Contact residues", fontsize=9)

        # legend patches
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor="#2e7d32", label="< 3 Å (close)"),
            Patch(facecolor="#f9a825", label="3–4 Å"),
            Patch(facecolor="#e64a19", label="> 4 Å"),
        ]
        ax_bar.legend(handles=legend_elements, fontsize=7, loc="lower right")
    else:
        ax_bar.text(0.5, 0.5, "No contact residues provided",
                    ha="center", va="center", fontsize=10, transform=ax_bar.transAxes)
        ax_bar.axis("off")

    fig.suptitle(
        f"Interaction map — {Path(docked_sdf).stem} (pose {pose_index})",
        fontsize=10, y=1.01,
    )
    fig.tight_layout()
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output_png


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
