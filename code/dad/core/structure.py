"""
dad.core.structure — Stage 4: Protein structure prediction via ColabFold / AlphaFold2 / AlphaFold3.

Wraps ColabFold v1.6.1 (AF2 backend, default) and parses AF3 CIF output.
All default parameters are transcribed from AlphaFold2.ipynb seed notebook.

Key defaults (from seed notebook)
----------------------------------
model_type        = "auto"               (ptm for monomer, multimer_v3 for complex)
msa_mode          = "mmseqs2_uniref_env"
num_recycles      = 3
num_seeds         = 1
num_models        = 5
use_dropout       = False
use_amber         = False
template_mode     = "none"
pair_mode         = "unpaired_paired"
pairing_strategy  = "greedy"
plddt_flag_thresh = 70.0

AF3 support
-----------
AlphaFold3 (Google DeepMind web server) produces multi-model CIF files.
Use ``load_af3_results()`` to parse an AF3 job output directory into
``StructurePrediction`` objects without re-running ColabFold.

Implementation owner: Mr_Struct (Phase B).
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional

from dad.io import ProteinInput, TriageDecision, StructurePrediction


# ─────────────────────────────────────────────────────────────────────────────
# Main prediction entry point
# ─────────────────────────────────────────────────────────────────────────────

def predict_structures(
    proteins: List[ProteinInput],
    decisions: List[TriageDecision],
    config: Dict,
    output_dir: str = "results/structures",
) -> List[StructurePrediction]:
    """Predict 3D structures for all triage-accepted ORFs.

    Only proteins with ``TriageDecision.verdict in {"accept", "downrank"}``
    are processed. For proteins with a dock_region defined, only that
    subsequence is submitted to ColabFold (reduces compute and avoids
    low-pLDDT TM helices polluting the docking region).

    Parameters
    ----------
    proteins : list of ProteinInput
        Full list of ORF translation records.
    decisions : list of TriageDecision
        Triage verdicts; used to filter and trim sequences.
    config : dict
        DAD config dict (see ``config.yaml``).
    output_dir : str, optional
        Base directory for all ColabFold job output folders.

    Returns
    -------
    list of StructurePrediction
        One result per successfully predicted structure.
        Structures with ``mean_plddt < config['structure']['plddt_min']``
        have ``low_confidence_flag = True`` but are still returned.

    Raises
    ------
    RuntimeError
        If ColabFold subprocess fails for a protein.
    """
    struct_cfg = config.get("structure", {})
    plddt_min = struct_cfg.get("plddt_min", 70.0)

    verdict_map: Dict[str, TriageDecision] = {d.seq_id: d for d in decisions}

    results: List[StructurePrediction] = []
    for prot in proteins:
        dec = verdict_map.get(prot.seq_id)
        if dec is None or dec.verdict == "exclude":
            continue

        seq = prot.sequence
        if dec.dock_region_start is not None and dec.dock_region_end is not None:
            seq = seq[dec.dock_region_start - 1: dec.dock_region_end]

        job_dir = str(Path(output_dir) / prot.seq_id)
        cf_cfg = struct_cfg
        cf_out = run_colabfold(
            seq_id=prot.seq_id,
            sequence=seq,
            job_dir=job_dir,
            model_type=cf_cfg.get("model_type", "auto"),
            msa_mode=cf_cfg.get("msa_mode", "mmseqs2_uniref_env"),
            num_recycles=cf_cfg.get("num_recycles", 3),
            num_seeds=cf_cfg.get("num_seeds", 1),
            num_models=cf_cfg.get("num_models", 5),
            use_dropout=cf_cfg.get("use_dropout", False),
            use_amber=cf_cfg.get("use_amber", False),
            template_mode=cf_cfg.get("template_mode", "none"),
            pair_mode=cf_cfg.get("pair_mode", "unpaired_paired"),
            pairing_strategy=cf_cfg.get("pairing_strategy", "greedy"),
            save_all=cf_cfg.get("save_all", False),
            save_recycles=cf_cfg.get("save_recycles", False),
        )

        rank1_pdb = _find_rank1_pdb(job_dir, prot.seq_id)
        quality = assess_structure_quality(rank1_pdb, plddt_min=plddt_min)

        results.append(StructurePrediction(
            seq_id=prot.seq_id,
            pdb_path=rank1_pdb,
            mean_plddt=quality["mean_plddt"],
            plddt_per_residue=quality["plddt_per_residue"],
            model_type=cf_out.get("model_type"),
            num_recycles=cf_cfg.get("num_recycles", 3),
            low_confidence_flag=quality["low_confidence_flag"],
            job_dir=job_dir,
        ))

    return results


def run_colabfold(
    seq_id: str,
    sequence: str,
    job_dir: str,
    model_type: str = "auto",
    msa_mode: str = "mmseqs2_uniref_env",
    num_recycles: int = 3,
    num_seeds: int = 1,
    num_models: int = 5,
    use_dropout: bool = False,
    use_amber: bool = False,
    template_mode: str = "none",
    pair_mode: str = "unpaired_paired",
    pairing_strategy: str = "greedy",
    save_all: bool = False,
    save_recycles: bool = False,
    user_agent: str = "dad/colabfold",
) -> Dict:
    """Run a single ColabFold prediction job.

    This is a thin wrapper around the ``colabfold.batch.run`` API,
    replicating the exact call sequence from ``AlphaFold2.ipynb``.

    Lazy-imports ``colabfold`` so the module can be imported on CPU-only
    machines (e.g., during Snakemake dry-run) without error.

    Parameters
    ----------
    seq_id : str
        Job name / identifier.
    sequence : str
        Amino acid sequence (use ``:`` separator for complexes).
    job_dir : str
        Directory where ColabFold writes all result files.
    model_type : str, optional
        ColabFold model type. ``"auto"`` selects ptm/multimer_v3 based on input.
    msa_mode : str, optional
        MSA generation mode. ``"mmseqs2_uniref_env"`` is the default.
    num_recycles : int, optional
        Number of AlphaFold2 recycling steps (default 3).
    num_seeds : int, optional
        Number of random seeds for stochastic sampling.
    num_models : int, optional
        Number of AF2 models to run (1–5).
    use_dropout : bool, optional
        Enable MC dropout for uncertainty sampling.
    use_amber : bool, optional
        Run AMBER relaxation on top-ranked structure.
    template_mode : str, optional
        Template usage: ``"none"``, ``"pdb100"``, or ``"custom"``.
    pair_mode : str, optional
        MSA pairing mode for complexes.
    pairing_strategy : str, optional
        Pairing strategy: ``"greedy"`` or ``"complete"``.
    save_all : bool, optional
        Save all model PDBs (not just rank_1).
    save_recycles : bool, optional
        Save intermediate recycle structures.
    user_agent : str, optional
        User-agent string for ColabFold MMseqs2 API calls.

    Returns
    -------
    dict
        ColabFold ``results`` dict containing ``rank``, ``metric``, etc.

    Raises
    ------
    ImportError
        If ``colabfold`` package is not installed (Colab/GPU environment required).
    RuntimeError
        If ColabFold job fails.
    """
    # Lazy import — only available in ColabFold/GPU environment
    try:
        from colabfold.batch import get_queries, run, set_model_type
        from colabfold.download import download_alphafold_params
        from colabfold.utils import setup_logging
    except ImportError as exc:
        raise ImportError(
            "colabfold is not installed. Run this function inside a ColabFold "
            "Colab environment or install via: "
            "pip install colabfold[alphafold-minus-jax]"
        ) from exc

    import warnings
    import logging
    warnings.filterwarnings("ignore", category=FutureWarning)

    job_path = Path(job_dir)
    job_path.mkdir(parents=True, exist_ok=True)

    queries_path = job_path / f"{seq_id}.csv"
    queries_path.write_text(f"id,sequence\n{seq_id},{sequence}\n", encoding="utf-8")

    use_templates = template_mode in ("pdb100", "custom")
    custom_template_path = None

    log_path = job_path / "log.txt"
    setup_logging(log_path)

    queries, is_complex = get_queries(str(queries_path))
    resolved_model_type = set_model_type(is_complex, model_type)

    use_cluster_profile = "multimer" not in resolved_model_type

    num_relax = 1 if use_amber else 0

    download_alphafold_params(resolved_model_type, Path("."))

    results = run(
        queries=queries,
        result_dir=str(job_path),
        use_templates=use_templates,
        custom_template_path=custom_template_path,
        num_relax=num_relax,
        msa_mode=msa_mode,
        model_type=resolved_model_type,
        num_models=num_models,
        num_recycles=num_recycles,
        recycle_early_stop_tolerance=None,
        relax_max_iterations=200,
        num_seeds=num_seeds,
        use_dropout=use_dropout,
        model_order=[1, 2, 3, 4, 5],
        is_complex=is_complex,
        data_dir=Path("."),
        keep_existing_results=False,
        rank_by="auto",
        pair_mode=pair_mode,
        pairing_strategy=pairing_strategy,
        stop_at_score=100.0,
        prediction_callback=None,
        dpi=200,
        zip_results=False,
        save_all=save_all,
        max_msa=None,
        use_cluster_profile=use_cluster_profile,
        input_features_callback=None,
        save_recycles=save_recycles,
        user_agent=user_agent,
    )

    results["model_type"] = resolved_model_type
    return results


# ─────────────────────────────────────────────────────────────────────────────
# AF3 support — parse AlphaFold3 Server output without re-running ColabFold
# ─────────────────────────────────────────────────────────────────────────────

def load_af3_results(
    af3_job_dir: str,
    seq_id: str,
    plddt_min: float = 70.0,
    model_index: int = 0,
) -> StructurePrediction:
    """Load an AlphaFold3 job output directory into a StructurePrediction object.

    AF3 server produces multi-model CIF files named
    ``<job>_model_{0..4}.cif`` with confidence JSON at
    ``<job>_summary_confidences_{0..4}.json``.

    Parameters
    ----------
    af3_job_dir : str
        Path to the AF3 output directory (contains ``*_model_*.cif`` files).
    seq_id : str
        Protein identifier to attach to the result.
    plddt_min : float, optional
        pLDDT threshold for ``low_confidence_flag``.
    model_index : int, optional
        Which model to select (0 = highest-confidence model in AF3 server output).

    Returns
    -------
    StructurePrediction
        Populated from the AF3 CIF and summary JSON.

    Raises
    ------
    FileNotFoundError
        If no CIF or confidence JSON files are found in ``af3_job_dir``.
    """
    job_path = Path(af3_job_dir)

    cif_files = sorted(job_path.glob(f"*_model_{model_index}.cif"))
    if not cif_files:
        cif_files = sorted(job_path.glob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files found in {af3_job_dir}")
    cif_path = cif_files[0]

    conf_files = sorted(job_path.glob(f"*_summary_confidences_{model_index}.json"))
    if not conf_files:
        conf_files = sorted(job_path.glob("*_summary_confidences_*.json"))

    mean_plddt = 0.0
    plddt_per_residue: List[float] = []
    ptm_score = None

    if conf_files:
        with open(conf_files[0], encoding="utf-8") as fh:
            conf = json.load(fh)
        atom_plddt = conf.get("atom_plddts") or conf.get("plddt") or []
        if atom_plddt:
            plddt_per_residue = _atom_plddt_to_residue(atom_plddt, cif_path)
            mean_plddt = sum(plddt_per_residue) / len(plddt_per_residue) if plddt_per_residue else 0.0
        ptm_score = conf.get("ptm") or conf.get("iptm")

    pdb_path = _convert_cif_to_pdb(str(cif_path))

    return StructurePrediction(
        seq_id=seq_id,
        pdb_path=pdb_path,
        mean_plddt=mean_plddt,
        ptm_score=ptm_score,
        plddt_per_residue=plddt_per_residue,
        model_type="alphafold3",
        num_recycles=None,
        low_confidence_flag=mean_plddt < plddt_min,
        job_dir=af3_job_dir,
    )


def load_af3_all_models(
    af3_job_dir: str,
    seq_id: str,
    plddt_min: float = 70.0,
) -> List[StructurePrediction]:
    """Load all 5 AF3 models from a job directory.

    Parameters
    ----------
    af3_job_dir : str
        Path to the AF3 output directory.
    seq_id : str
        Protein identifier.
    plddt_min : float, optional
        pLDDT threshold for ``low_confidence_flag``.

    Returns
    -------
    list of StructurePrediction
        One per model (model_0 … model_4), sorted by mean pLDDT descending.
    """
    results = []
    job_path = Path(af3_job_dir)
    for model_idx in range(5):
        cif_files = sorted(job_path.glob(f"*_model_{model_idx}.cif"))
        if not cif_files:
            continue
        try:
            sp = load_af3_results(af3_job_dir, seq_id, plddt_min, model_index=model_idx)
            results.append(sp)
        except FileNotFoundError:
            continue
    results.sort(key=lambda x: x.mean_plddt, reverse=True)
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Structure quality assessment
# ─────────────────────────────────────────────────────────────────────────────

def assess_structure_quality(
    pdb_path: str,
    plddt_min: float = 70.0,
) -> Dict:
    """Parse pLDDT scores from a ColabFold PDB file and assess quality.

    pLDDT values are stored in the B-factor column of ColabFold PDB output.

    Parameters
    ----------
    pdb_path : str
        Path to the ColabFold PDB file.
    plddt_min : float, optional
        Threshold below which ``low_confidence_flag`` is set to True.

    Returns
    -------
    dict
        Keys: ``mean_plddt`` (float), ``plddt_per_residue`` (list of float),
        ``low_confidence_flag`` (bool), ``high_confidence_fraction`` (float).

    Raises
    ------
    FileNotFoundError
        If ``pdb_path`` does not exist.
    """
    pdb_file = Path(pdb_path)
    if not pdb_file.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    residue_plddt: Dict[tuple, float] = {}
    with open(pdb_file, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                chain = line[21]
                res_seq = int(line[22:26].strip())
                b_factor = float(line[60:66].strip())
            except (ValueError, IndexError):
                continue
            key = (chain, res_seq)
            # keep first atom's B-factor (CA or first ATOM) per residue
            if key not in residue_plddt:
                residue_plddt[key] = b_factor

    if not residue_plddt:
        return {
            "mean_plddt": 0.0,
            "plddt_per_residue": [],
            "low_confidence_flag": True,
            "high_confidence_fraction": 0.0,
        }

    plddt_values = list(residue_plddt.values())
    mean_plddt = sum(plddt_values) / len(plddt_values)
    high_conf_frac = sum(1 for v in plddt_values if v >= plddt_min) / len(plddt_values)

    return {
        "mean_plddt": round(mean_plddt, 3),
        "plddt_per_residue": [round(v, 3) for v in plddt_values],
        "low_confidence_flag": mean_plddt < plddt_min,
        "high_confidence_fraction": round(high_conf_frac, 3),
    }


def select_dock_region_pdb(
    pdb_path: str,
    dock_region_start: Optional[int],
    dock_region_end: Optional[int],
    output_path: str,
) -> str:
    """Extract a sub-chain PDB for the dock-competent region.

    For membrane proteins where only the extracellular domain is docked,
    this function trims the PDB to the specified residue range.
    If ``dock_region_start`` and ``dock_region_end`` are both None,
    the full PDB is copied to ``output_path`` unchanged.

    Parameters
    ----------
    pdb_path : str
        Input PDB file path (full structure).
    dock_region_start : int or None
        First residue index (1-based) to retain.
    dock_region_end : int or None
        Last residue index (1-based) to retain.
    output_path : str
        Destination path for the trimmed PDB.

    Returns
    -------
    str
        ``output_path`` on success.

    Raises
    ------
    FileNotFoundError
        If ``pdb_path`` does not exist.
    """
    src = Path(pdb_path)
    if not src.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    if dock_region_start is None or dock_region_end is None:
        shutil.copy2(str(src), output_path)
        return output_path

    # Use Bio.PDB if available, otherwise pure-text fallback
    try:
        from Bio.PDB import PDBParser, PDBIO, Select

        class _RegionSelect(Select):
            def accept_residue(self, residue):
                res_id = residue.get_id()[1]
                return dock_region_start <= res_id <= dock_region_end

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("target", str(src))
        io = PDBIO()
        io.set_structure(structure)
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        io.save(output_path, _RegionSelect())

    except ImportError:
        # Fallback: plain-text residue range filter
        out_lines = []
        with open(src, encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(("ATOM  ", "HETATM")):
                    try:
                        res_seq = int(line[22:26].strip())
                    except ValueError:
                        out_lines.append(line)
                        continue
                    if dock_region_start <= res_seq <= dock_region_end:
                        out_lines.append(line)
                else:
                    out_lines.append(line)
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as fh:
            fh.writelines(out_lines)

    return output_path


# ─────────────────────────────────────────────────────────────────────────────
# AW1_ref reuse helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_aw1_asset(asset_id: str, manifest_path: Path) -> Path:
    """Resolve a canonical AW1_ref file path via manifest, not filename inference.

    Reads ``06_Report/Mr_Repro/aw1_ref_manifest.tsv`` (or any TSV sharing the
    same column layout) and returns the ``path`` field for the given
    ``asset_id``.

    Manifest columns (tab-separated, header on row 1):
        asset_id  seq_id  gene  stage  asset_type  path  selection_reason  sha256

    Parameters
    ----------
    asset_id : str
        The manifest ``asset_id`` value (e.g. ``"MCP_AF2_PDB"``).
    manifest_path : Path
        Path to the manifest TSV file.

    Returns
    -------
    Path
        Resolved absolute path to the asset file.

    Raises
    ------
    FileNotFoundError
        If ``manifest_path`` does not exist.
    KeyError
        If ``asset_id`` is not found in the manifest.
    FileNotFoundError
        If the resolved path does not exist on disk.
    """
    if not manifest_path.exists():
        raise FileNotFoundError(f"AW1_ref manifest not found: {manifest_path}")

    import csv as _csv
    with open(manifest_path, newline="", encoding="utf-8-sig") as fh:
        reader = _csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("asset_id", "").strip() == asset_id:
                raw_path = row.get("path", "").strip()
                resolved = Path(raw_path)
                if not resolved.is_absolute():
                    # Paths in manifest are relative to DAD project root
                    resolved = manifest_path.parent.parent.parent / raw_path
                resolved = resolved.resolve()
                if not resolved.exists():
                    raise FileNotFoundError(
                        f"Manifest asset '{asset_id}' path does not exist: {resolved}"
                    )
                return resolved

    raise KeyError(f"asset_id '{asset_id}' not found in manifest: {manifest_path}")


def load_existing_pdb(
    pdb_path: str,
    seq_id: str,
    plddt_min: float = 70.0,
    manifest_path: Optional[Path] = None,
    asset_id: Optional[str] = None,
) -> StructurePrediction:
    """Wrap an existing PDB (e.g., from AW1_ref) into a StructurePrediction.

    Use this to reuse Stage 4 outputs from AW1_ref without re-running ColabFold.

    When ``manifest_path`` and ``asset_id`` are both provided, the canonical
    path is resolved via ``load_aw1_asset()`` first.  The ``pdb_path`` argument
    is used as a fallback if the manifest is absent or the asset_id is not found.

    Parameters
    ----------
    pdb_path : str
        Fallback path to an existing PDB file (ColabFold or AF3-converted).
    seq_id : str
        Protein identifier.
    plddt_min : float, optional
        pLDDT threshold for ``low_confidence_flag``.
    manifest_path : Path, optional
        Path to ``aw1_ref_manifest.tsv``.  If provided together with
        ``asset_id``, the manifest entry takes priority over ``pdb_path``.
    asset_id : str, optional
        Manifest ``asset_id`` key (e.g. ``"MCP_AF2_PDB"``).

    Returns
    -------
    StructurePrediction
    """
    resolved_path = pdb_path
    if manifest_path is not None and asset_id is not None:
        try:
            resolved_path = str(load_aw1_asset(asset_id, manifest_path))
        except (FileNotFoundError, KeyError):
            # Manifest not yet available or asset not listed — fall back to direct path
            resolved_path = pdb_path

    quality = assess_structure_quality(resolved_path, plddt_min=plddt_min)
    return StructurePrediction(
        seq_id=seq_id,
        pdb_path=str(Path(resolved_path).resolve()),
        mean_plddt=quality["mean_plddt"],
        plddt_per_residue=quality["plddt_per_residue"],
        model_type="precomputed",
        low_confidence_flag=quality["low_confidence_flag"],
        job_dir=str(Path(resolved_path).parent),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _find_rank1_pdb(job_dir: str, seq_id: str) -> str:
    """Return the path to the rank_1 PDB produced by ColabFold."""
    job_path = Path(job_dir)
    # ColabFold v1.6.1 naming: <jobname>_unrelaxed_rank_001_*.pdb
    candidates = sorted(job_path.glob("*_unrelaxed_rank_001_*.pdb"))
    if candidates:
        return str(candidates[0])
    # Legacy naming: *_rank_1_*.pdb
    candidates = sorted(job_path.glob("*_rank_1_*.pdb"))
    if candidates:
        return str(candidates[0])
    # Fallback: first PDB in directory
    candidates = sorted(job_path.glob("*.pdb"))
    if candidates:
        return str(candidates[0])
    raise FileNotFoundError(f"No rank_1 PDB found in {job_dir}")


def _convert_cif_to_pdb(cif_path: str) -> str:
    """Convert an mmCIF file to PDB format, returning the PDB path.

    Uses Bio.PDB if available, otherwise attempts gemmi CLI, then keeps CIF
    and notes the path (downstream code must handle CIF-native reading).
    """
    out_path = cif_path.replace(".cif", ".pdb")
    if Path(out_path).exists():
        return out_path

    try:
        from Bio.PDB import MMCIFParser, PDBIO
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("af3", cif_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(out_path)
        return out_path
    except ImportError:
        pass

    # gemmi CLI fallback (best-effort — not required)
    try:
        ret = subprocess.run(
            ["gemmi", "convert", "--to=pdb", cif_path, out_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if ret.returncode == 0 and Path(out_path).exists():
            return out_path
    except (FileNotFoundError, OSError):
        pass

    # Return cif path as last resort — caller must handle
    return cif_path


def _atom_plddt_to_residue(atom_plddt: List[float], cif_path: Path) -> List[float]:
    """Reduce per-atom pLDDT array to per-residue by taking the CA (or first atom).

    AF3 ``summary_confidences`` stores atom-level pLDDT.
    We collapse to one value per residue using the first atom encountered.
    """
    residue_plddt: List[float] = []
    prev_res = None
    idx = 0

    try:
        with open(cif_path, encoding="utf-8", errors="ignore") as fh:
            in_atom_loop = False
            col_map: Dict[str, int] = {}
            col_idx = 0
            for line in fh:
                line = line.strip()
                if line == "_atom_site.label_seq_id":
                    in_atom_loop = True
                    col_map["seq_id"] = col_idx
                elif line.startswith("_atom_site."):
                    if in_atom_loop:
                        col_idx += 1
                elif in_atom_loop and line and not line.startswith("_") and not line.startswith("#"):
                    parts = line.split()
                    seq_col = col_map.get("seq_id", 6)
                    try:
                        res_id = int(parts[seq_col])
                    except (IndexError, ValueError):
                        continue
                    if res_id != prev_res:
                        if idx < len(atom_plddt):
                            residue_plddt.append(atom_plddt[idx])
                        prev_res = res_id
                    idx += 1
    except Exception:
        pass

    if not residue_plddt and atom_plddt:
        # Simple fallback: every 7 atoms per residue (rough average)
        residue_plddt = atom_plddt[::7]

    return residue_plddt
