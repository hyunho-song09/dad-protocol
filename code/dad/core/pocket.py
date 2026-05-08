"""
dad.core.pocket — Stage 5: Binding pocket detection via P2Rank.

Wraps P2Rank with AlphaFold-specific configuration (``-c alphafold``).
Returns ranked pocket centroids for use as GNINA docking box centers.

AW1_ref reuse
-------------
``load_existing_predictions_csv()`` parses an existing P2Rank
``*_predictions.csv`` (e.g., ``AW1_ref/structure.pdb_predictions_MCP.csv``)
directly into ``PocketResult`` objects without re-running P2Rank.

CSV column layout (space-padded, comma-delimited):
  name, rank, score, probability, sas_points, surf_atoms,
  center_x, center_y, center_z, residue_ids, surf_atom_ids

Key defaults (from DAD_project_plan.md Stage 5)
------------------------------------------------
tool           = "p2rank"
alphafold_mode = True   (adds -c alphafold flag)
top_n_pockets  = 3
min_score      = 0.0

Implementation owner: Mr_Struct (Phase B).
"""

from __future__ import annotations

import csv
import io
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from dad.core.structure import load_aw1_asset  # manifest resolver

from dad.io import StructurePrediction, PocketResult


# ─────────────────────────────────────────────────────────────────────────────
# Main pocket detection entry point
# ─────────────────────────────────────────────────────────────────────────────

def detect_pockets(
    structures: List[StructurePrediction],
    config: Dict,
    output_dir: str = "results/pockets",
) -> Dict[str, List[PocketResult]]:
    """Detect binding pockets for all predicted structures using P2Rank.

    Parameters
    ----------
    structures : list of StructurePrediction
        Structure prediction results from Stage 4.
    config : dict
        DAD config dict (see ``config.yaml``). Reads ``config['pocket']``.
    output_dir : str, optional
        Base directory for P2Rank output files.

    Returns
    -------
    dict
        Mapping ``seq_id -> list of PocketResult`` sorted by P2Rank score
        (descending). Only the top-N pockets (``config['pocket']['top_n_pockets']``)
        per protein are returned.

    Raises
    ------
    RuntimeError
        If P2Rank subprocess fails.
    """
    pocket_cfg = config.get("pocket", {})
    p2rank_path = pocket_cfg.get("p2rank_path", "prank")
    alphafold_mode = pocket_cfg.get("alphafold_mode", True)
    top_n = pocket_cfg.get("top_n_pockets", 3)
    min_score = pocket_cfg.get("min_score", 0.0)

    out_base = Path(output_dir)
    results: Dict[str, List[PocketResult]] = {}

    for struct in structures:
        job_out = str(out_base / struct.seq_id)
        try:
            predictions_csv = run_p2rank(
                pdb_path=struct.pdb_path,
                output_dir=job_out,
                p2rank_path=p2rank_path,
                alphafold_mode=alphafold_mode,
                top_n=top_n,
            )
            pockets = parse_p2rank_output(
                predictions_csv=predictions_csv,
                seq_id=struct.seq_id,
                top_n=top_n,
                min_score=min_score,
            )
        except (FileNotFoundError, RuntimeError) as exc:
            raise RuntimeError(
                f"P2Rank failed for {struct.seq_id}: {exc}"
            ) from exc

        results[struct.seq_id] = pockets

    return results


def run_p2rank(
    pdb_path: str,
    output_dir: str,
    p2rank_path: str = "prank",
    alphafold_mode: bool = True,
    top_n: int = 3,
) -> str:
    """Run P2Rank on a single PDB file.

    Parameters
    ----------
    pdb_path : str
        Absolute path to the input PDB file.
    output_dir : str
        Directory where P2Rank writes its output.
    p2rank_path : str, optional
        Path to the ``prank`` executable.
    alphafold_mode : bool, optional
        If True, adds ``-c alphafold`` to the P2Rank command.
        This uses AF2-optimized pocket scoring parameters.
    top_n : int, optional
        Number of top pockets to return (applied at parsing stage).

    Returns
    -------
    str
        Path to the P2Rank ``predictions.csv`` output file.

    Raises
    ------
    FileNotFoundError
        If ``pdb_path`` is not found.
    RuntimeError
        If P2Rank subprocess returns non-zero exit code.
    """
    pdb_file = Path(pdb_path)
    if not pdb_file.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    cmd = [p2rank_path, "predict", "-f", str(pdb_file), "-o", str(out_path)]
    if alphafold_mode:
        cmd += ["-c", "alphafold"]

    ret = subprocess.run(cmd, capture_output=True, text=True)
    if ret.returncode != 0:
        raise RuntimeError(
            f"P2Rank failed (exit {ret.returncode}):\n{ret.stderr}"
        )

    # P2Rank writes <pdb_stem>_predictions.csv
    pdb_stem = pdb_file.stem
    csv_candidates = [
        out_path / f"{pdb_stem}_predictions.csv",
        out_path / "predictions.csv",
        out_path / f"{pdb_file.name}_predictions.csv",
    ]
    for candidate in csv_candidates:
        if candidate.exists():
            return str(candidate)

    # Glob fallback
    hits = list(out_path.glob("*predictions*.csv"))
    if hits:
        return str(hits[0])

    raise FileNotFoundError(
        f"P2Rank ran but predictions.csv not found in {output_dir}. "
        f"P2Rank stdout: {ret.stdout}"
    )


def parse_p2rank_output(
    predictions_csv: str,
    seq_id: str,
    top_n: int = 3,
    min_score: float = 0.0,
) -> List[PocketResult]:
    """Parse P2Rank ``predictions.csv`` into a list of PocketResult objects.

    Handles both AW1_ref format (space-padded, comma-delimited) and standard
    P2Rank CSV output.

    P2Rank output CSV columns (space-padded, comma-separated):
    ``name, rank, score, probability, sas_points, surf_atoms,
    center_x, center_y, center_z, residue_ids, surf_atom_ids``

    Parameters
    ----------
    predictions_csv : str
        Path to the P2Rank ``predictions.csv`` (or ``structure.pdb_predictions_*.csv``).
    seq_id : str
        Protein identifier to attach to each PocketResult.
    top_n : int, optional
        Maximum number of pockets to return.
    min_score : float, optional
        Minimum pocket score threshold (pockets below are discarded).

    Returns
    -------
    list of PocketResult
        Parsed pocket records sorted by rank. Length <= ``top_n``.

    Raises
    ------
    FileNotFoundError
        If ``predictions_csv`` does not exist.
    ValueError
        If CSV cannot be parsed (unexpected format).
    """
    csv_file = Path(predictions_csv)
    if not csv_file.exists():
        raise FileNotFoundError(f"P2Rank predictions CSV not found: {predictions_csv}")

    raw = csv_file.read_text(encoding="utf-8")

    # Strip leading/trailing whitespace from each field (P2Rank pads with spaces)
    cleaned_lines = []
    for line in raw.splitlines():
        if not line.strip() or line.strip().startswith("#"):
            continue
        cleaned_lines.append(line)

    if len(cleaned_lines) < 2:
        return []

    reader = csv.DictReader(
        io.StringIO("\n".join(cleaned_lines)),
        skipinitialspace=True,
    )
    # Normalise header names (strip whitespace)
    fieldnames = [f.strip() for f in (reader.fieldnames or [])]
    reader.fieldnames = fieldnames

    pockets: List[PocketResult] = []
    for row in reader:
        row = {k.strip(): v.strip() for k, v in row.items() if k}
        try:
            rank = int(row.get("rank", 0))
            score = float(row.get("score", 0.0))
            center_x = float(row.get("center_x", 0.0))
            center_y = float(row.get("center_y", 0.0))
            center_z = float(row.get("center_z", 0.0))
            surf_atoms = int(row.get("surf_atoms", 0)) if row.get("surf_atoms") else None
        except (ValueError, KeyError) as exc:
            raise ValueError(
                f"Cannot parse P2Rank CSV row {row}: {exc}"
            ) from exc

        if score < min_score:
            continue

        residues = _parse_residue_ids(row.get("residue_ids", ""))

        pockets.append(PocketResult(
            seq_id=seq_id,
            pocket_rank=rank,
            score=score,
            center_x=center_x,
            center_y=center_y,
            center_z=center_z,
            residues=residues,
            surf_atoms=surf_atoms,
        ))

    pockets.sort(key=lambda p: p.pocket_rank)
    return pockets[:top_n]


# ─────────────────────────────────────────────────────────────────────────────
# AW1_ref reuse — load existing predictions CSV directly
# ─────────────────────────────────────────────────────────────────────────────

def load_existing_predictions_csv(
    predictions_csv: str,
    seq_id: str,
    top_n: int = 3,
    min_score: float = 0.0,
    manifest_path: Optional[Path] = None,
    asset_id: Optional[str] = None,
) -> List[PocketResult]:
    """Load an existing P2Rank predictions CSV from AW1_ref or prior runs.

    This is the primary path for AW1_ref reuse (Phase B principle:
    reuse existing Stage 5 outputs before re-running P2Rank).

    When ``manifest_path`` and ``asset_id`` are both provided, the canonical
    path is resolved via ``load_aw1_asset()`` first.  The ``predictions_csv``
    argument is used as a fallback if the manifest is absent or the asset is
    not found.

    Accepts both ``structure.pdb_predictions_MCP.csv`` naming
    (AW1_ref top-level) and ``structure.cif_predictions.csv`` naming
    (per-protein P2Rank web output).

    Parameters
    ----------
    predictions_csv : str
        Fallback path to the existing predictions CSV file.
    seq_id : str
        Protein identifier to assign to parsed pockets.
    top_n : int, optional
        Maximum number of pockets to return.
    min_score : float, optional
        Minimum pocket score threshold.
    manifest_path : Path, optional
        Path to ``aw1_ref_manifest.tsv``.  If provided together with
        ``asset_id``, the manifest entry takes priority over ``predictions_csv``.
    asset_id : str, optional
        Manifest ``asset_id`` key (e.g. ``"MCP_P2RANK"``).

    Returns
    -------
    list of PocketResult
    """
    resolved_csv = predictions_csv
    if manifest_path is not None and asset_id is not None:
        try:
            resolved_csv = str(load_aw1_asset(asset_id, manifest_path))
        except (FileNotFoundError, KeyError):
            resolved_csv = predictions_csv

    return parse_p2rank_output(resolved_csv, seq_id, top_n=top_n, min_score=min_score)


def pocket_to_box(
    pocket: PocketResult,
    ligand_max_dim: float,
    box_size_min: float = 22.0,
    box_padding: float = 10.0,
) -> Dict:
    """Compute GNINA docking box dimensions from a pocket centroid.

    Box size rule (DAD Stage 7, zero-tuning):
        side = max(box_size_min, ligand_max_dim + box_padding)

    All three axes use the same side length for simplicity.

    Parameters
    ----------
    pocket : PocketResult
        P2Rank pocket result providing the box center coordinates.
    ligand_max_dim : float
        Maximum atom-pair distance in the prepared ligand conformer (A).
    box_size_min : float, optional
        Minimum box side length in A (default 22, from seed notebook).
    box_padding : float, optional
        Padding added to ``ligand_max_dim`` (default 10 A).

    Returns
    -------
    dict
        Keys: ``center_x``, ``center_y``, ``center_z``,
        ``size_x``, ``size_y``, ``size_z`` (all floats).

    Raises
    ------
    ValueError
        If ``ligand_max_dim`` is negative.
    """
    if ligand_max_dim < 0:
        raise ValueError(f"ligand_max_dim must be non-negative, got {ligand_max_dim}")
    side = max(box_size_min, ligand_max_dim + box_padding)
    return {
        "center_x": pocket.center_x,
        "center_y": pocket.center_y,
        "center_z": pocket.center_z,
        "size_x": side,
        "size_y": side,
        "size_z": side,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _parse_residue_ids(residue_ids_str: str) -> List[int]:
    """Parse P2Rank residue_ids field into a list of integer residue numbers.

    Input examples:
      "A_122 A_126 A_129"  -> [122, 126, 129]
      "122 126 129"        -> [122, 126, 129]
    """
    if not residue_ids_str:
        return []
    residues = []
    for token in residue_ids_str.split():
        token = token.strip()
        if "_" in token:
            token = token.split("_", 1)[-1]
        try:
            residues.append(int(token))
        except ValueError:
            continue
    return residues
