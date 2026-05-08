"""
Tests for dad.core.structure and dad.core.pocket.

All tests run on CPU without ColabFold or P2Rank installed.
They exercise:
  - assess_structure_quality() against AW1_ref/structure_MCP.pdb
  - load_existing_pdb() AW1_ref reuse path
  - load_af3_results() against AW1_ref AF3 CIF output
  - select_dock_region_pdb() trimming
  - parse_p2rank_output() against AW1_ref/structure.pdb_predictions_MCP.csv
  - load_existing_predictions_csv() convenience wrapper
  - pocket_to_box() box sizing logic
  - _parse_residue_ids() internal helper
"""

import sys
import os
import tempfile
from pathlib import Path

# Allow running from repo root or tests/ directly
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent / "06_Report" / "Mr_Pipeline"))

import pytest

from dad.core.structure import (
    assess_structure_quality,
    load_existing_pdb,
    load_aw1_asset,
    select_dock_region_pdb,
    load_af3_results,
    load_af3_all_models,
    _find_rank1_pdb,
)
from dad.core.pocket import (
    parse_p2rank_output,
    load_existing_predictions_csv,
    pocket_to_box,
    _parse_residue_ids,
)
from dad.io import StructurePrediction, PocketResult

# ─────────────────────────────────────────────────────────────────────────────
# Paths to AW1_ref assets
# ─────────────────────────────────────────────────────────────────────────────

AW1_REF = Path(__file__).resolve().parent.parent.parent.parent / "AW1_ref"
MCP_PDB = AW1_REF / "structure_MCP.pdb"
CRP_PDB = AW1_REF / "structure_Crp.pdb"
RBS_PDB = AW1_REF / "structure_Rbs.pdb"
MCP_PRED_CSV = AW1_REF / "structure.pdb_predictions_MCP.csv"
CRP_PRED_CSV = AW1_REF / "structure.pdb_predictions_Crp.csv"
RBS_PRED_CSV = AW1_REF / "structure.pdb_predictions_Rbs.csv"

AF3_MCP_DIR = AW1_REF / "MCP" / "NA23_RS02135" / "last" / \
    "fold_2025_12_27_12_50_na23_rs02135_sensing_domain"

ASSETS_AVAILABLE = MCP_PDB.exists()


# ─────────────────────────────────────────────────────────────────────────────
# assess_structure_quality
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_assess_structure_quality_mcp_returns_dict():
    result = assess_structure_quality(str(MCP_PDB))
    assert isinstance(result, dict)
    for key in ("mean_plddt", "plddt_per_residue", "low_confidence_flag", "high_confidence_fraction"):
        assert key in result, f"Missing key: {key}"


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_assess_structure_quality_mcp_plddt_range():
    result = assess_structure_quality(str(MCP_PDB))
    assert 0.0 <= result["mean_plddt"] <= 100.0, "mean_plddt out of [0, 100]"
    for v in result["plddt_per_residue"]:
        assert 0.0 <= v <= 100.0, f"Residue pLDDT {v} out of range"


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_assess_structure_quality_mcp_residues_nonzero():
    result = assess_structure_quality(str(MCP_PDB))
    assert len(result["plddt_per_residue"]) > 0, "No residues parsed"


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_assess_structure_quality_mcp_high_confidence():
    # MCP structure from AW1_ref should have decent quality (expected > 50%)
    result = assess_structure_quality(str(MCP_PDB), plddt_min=50.0)
    assert result["high_confidence_fraction"] >= 0.0


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_assess_structure_quality_crp():
    result = assess_structure_quality(str(CRP_PDB))
    assert result["mean_plddt"] > 0.0


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_assess_structure_quality_rbs():
    result = assess_structure_quality(str(RBS_PDB))
    assert result["mean_plddt"] > 0.0


def test_assess_structure_quality_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        assess_structure_quality("/nonexistent/path/to.pdb")


# ─────────────────────────────────────────────────────────────────────────────
# load_existing_pdb (AW1_ref reuse path)
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_load_existing_pdb_returns_structure_prediction():
    sp = load_existing_pdb(str(MCP_PDB), seq_id="MCP")
    assert isinstance(sp, StructurePrediction)
    assert sp.seq_id == "MCP"
    assert sp.model_type == "precomputed"
    assert Path(sp.pdb_path).exists()


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_load_existing_pdb_plddt_populated():
    sp = load_existing_pdb(str(MCP_PDB), seq_id="MCP")
    assert sp.mean_plddt > 0.0
    assert len(sp.plddt_per_residue) > 0


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_load_existing_pdb_low_confidence_flag_correct():
    sp = load_existing_pdb(str(MCP_PDB), seq_id="MCP", plddt_min=100.0)
    # With threshold 100, everything should be low-confidence
    assert sp.low_confidence_flag is True

    sp2 = load_existing_pdb(str(MCP_PDB), seq_id="MCP", plddt_min=0.0)
    assert sp2.low_confidence_flag is False


# ─────────────────────────────────────────────────────────────────────────────
# select_dock_region_pdb
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_select_dock_region_none_copies_unchanged(tmp_path):
    out = str(tmp_path / "full.pdb")
    result = select_dock_region_pdb(str(MCP_PDB), None, None, out)
    assert result == out
    assert Path(out).exists()
    # Sizes should be the same
    assert Path(out).stat().st_size == MCP_PDB.stat().st_size


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_select_dock_region_trim_reduces_residues(tmp_path):
    out = str(tmp_path / "trimmed.pdb")
    result = select_dock_region_pdb(str(MCP_PDB), 50, 150, out)
    assert result == out
    assert Path(out).exists()
    # Trimmed file should be smaller than original
    assert Path(out).stat().st_size < MCP_PDB.stat().st_size


def test_select_dock_region_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        select_dock_region_pdb("/no/such.pdb", 1, 100, str(tmp_path / "out.pdb"))


# ─────────────────────────────────────────────────────────────────────────────
# AF3 loading (if AF3 assets present)
# ─────────────────────────────────────────────────────────────────────────────

AF3_AVAILABLE = AF3_MCP_DIR.exists()


@pytest.mark.skipif(not AF3_AVAILABLE, reason="AF3 CIF assets not found")
def test_load_af3_results_returns_structure_prediction():
    sp = load_af3_results(str(AF3_MCP_DIR), seq_id="MCP_AF3")
    assert isinstance(sp, StructurePrediction)
    assert sp.seq_id == "MCP_AF3"
    assert sp.model_type == "alphafold3"
    # pdb_path may be .cif if Bio/gemmi not installed (fallback accepted)
    assert Path(sp.pdb_path).exists() or sp.pdb_path.endswith(".cif")


@pytest.mark.skipif(not AF3_AVAILABLE, reason="AF3 CIF assets not found")
def test_load_af3_all_models_returns_list():
    sps = load_af3_all_models(str(AF3_MCP_DIR), seq_id="MCP_AF3")
    assert len(sps) >= 1
    # Sorted by mean_plddt descending
    for i in range(len(sps) - 1):
        assert sps[i].mean_plddt >= sps[i + 1].mean_plddt


def test_load_af3_results_missing_dir_raises():
    with pytest.raises(FileNotFoundError):
        load_af3_results("/nonexistent/af3_dir", seq_id="X")


# ─────────────────────────────────────────────────────────────────────────────
# parse_p2rank_output / load_existing_predictions_csv
# ─────────────────────────────────────────────────────────────────────────────

PRED_CSV_AVAILABLE = MCP_PRED_CSV.exists()


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_parse_p2rank_output_mcp_returns_pockets():
    pockets = parse_p2rank_output(str(MCP_PRED_CSV), seq_id="MCP", top_n=3)
    assert len(pockets) > 0
    assert len(pockets) <= 3


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_parse_p2rank_output_pocket_result_types():
    pockets = parse_p2rank_output(str(MCP_PRED_CSV), seq_id="MCP", top_n=3)
    for p in pockets:
        assert isinstance(p, PocketResult)
        assert p.seq_id == "MCP"
        assert isinstance(p.pocket_rank, int)
        assert isinstance(p.score, float)
        assert isinstance(p.center_x, float)
        assert isinstance(p.center_y, float)
        assert isinstance(p.center_z, float)


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_parse_p2rank_output_mcp_sorted_by_rank():
    pockets = parse_p2rank_output(str(MCP_PRED_CSV), seq_id="MCP", top_n=5)
    ranks = [p.pocket_rank for p in pockets]
    assert ranks == sorted(ranks), "Pockets not sorted by rank"


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_parse_p2rank_output_top_n_respected():
    pockets = parse_p2rank_output(str(MCP_PRED_CSV), seq_id="MCP", top_n=2)
    assert len(pockets) <= 2


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_parse_p2rank_output_min_score_filter():
    # With a very high min_score, should return nothing
    pockets = parse_p2rank_output(str(MCP_PRED_CSV), seq_id="MCP", min_score=999.0)
    assert pockets == []


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_load_existing_predictions_csv_mcp():
    pockets = load_existing_predictions_csv(str(MCP_PRED_CSV), seq_id="MCP")
    assert len(pockets) > 0
    assert pockets[0].seq_id == "MCP"


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_parse_p2rank_output_crp():
    pockets = parse_p2rank_output(str(CRP_PRED_CSV), seq_id="Crp", top_n=3)
    assert len(pockets) > 0


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_parse_p2rank_output_rbs():
    pockets = parse_p2rank_output(str(RBS_PRED_CSV), seq_id="Rbs", top_n=3)
    assert len(pockets) > 0


def test_parse_p2rank_output_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        parse_p2rank_output("/nonexistent/predictions.csv", seq_id="X")


def test_parse_p2rank_output_empty_csv(tmp_path):
    empty_csv = tmp_path / "empty.csv"
    empty_csv.write_text("name, rank, score, probability, sas_points, surf_atoms, center_x, center_y, center_z, residue_ids, surf_atom_ids\n")
    pockets = parse_p2rank_output(str(empty_csv), seq_id="X")
    assert pockets == []


# ─────────────────────────────────────────────────────────────────────────────
# pocket_to_box
# ─────────────────────────────────────────────────────────────────────────────

def _make_pocket(cx=1.0, cy=2.0, cz=3.0, score=10.0, rank=1):
    return PocketResult(
        seq_id="test",
        pocket_rank=rank,
        score=score,
        center_x=cx,
        center_y=cy,
        center_z=cz,
    )


def test_pocket_to_box_minimum_size():
    pocket = _make_pocket()
    box = pocket_to_box(pocket, ligand_max_dim=5.0, box_size_min=22.0, box_padding=10.0)
    # side = max(22, 5+10) = 22
    assert box["size_x"] == 22.0
    assert box["size_y"] == 22.0
    assert box["size_z"] == 22.0


def test_pocket_to_box_padding_dominates():
    pocket = _make_pocket()
    box = pocket_to_box(pocket, ligand_max_dim=20.0, box_size_min=22.0, box_padding=10.0)
    # side = max(22, 20+10) = 30
    assert box["size_x"] == 30.0


def test_pocket_to_box_center_passthrough():
    pocket = _make_pocket(cx=5.5, cy=-3.2, cz=11.1)
    box = pocket_to_box(pocket, ligand_max_dim=5.0)
    assert box["center_x"] == 5.5
    assert box["center_y"] == -3.2
    assert box["center_z"] == 11.1


def test_pocket_to_box_zero_ligand_dim():
    pocket = _make_pocket()
    box = pocket_to_box(pocket, ligand_max_dim=0.0)
    assert box["size_x"] == 22.0


def test_pocket_to_box_negative_ligand_dim_raises():
    pocket = _make_pocket()
    with pytest.raises(ValueError):
        pocket_to_box(pocket, ligand_max_dim=-1.0)


# ─────────────────────────────────────────────────────────────────────────────
# _parse_residue_ids
# ─────────────────────────────────────────────────────────────────────────────

def test_parse_residue_ids_chain_underscore_format():
    result = _parse_residue_ids("A_122 A_126 A_129")
    assert result == [122, 126, 129]


def test_parse_residue_ids_plain_format():
    result = _parse_residue_ids("122 126 129")
    assert result == [122, 126, 129]


def test_parse_residue_ids_empty_string():
    result = _parse_residue_ids("")
    assert result == []


def test_parse_residue_ids_ignores_non_numeric():
    result = _parse_residue_ids("A_122 B_abc A_130")
    assert result == [122, 130]


# ─────────────────────────────────────────────────────────────────────────────
# Integration: load_existing_pdb + load_existing_predictions_csv + pocket_to_box
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.skipif(
    not (ASSETS_AVAILABLE and PRED_CSV_AVAILABLE),
    reason="AW1_ref assets not found",
)
def test_aw1ref_mcp_full_pipeline_no_gpu():
    """End-to-end AW1_ref reuse: load PDB -> load pockets -> compute box."""
    sp = load_existing_pdb(str(MCP_PDB), seq_id="MCP")
    assert sp.mean_plddt > 0.0

    pockets = load_existing_predictions_csv(str(MCP_PRED_CSV), seq_id="MCP", top_n=3)
    assert len(pockets) > 0

    top_pocket = pockets[0]
    box = pocket_to_box(top_pocket, ligand_max_dim=8.0)
    assert box["size_x"] >= 18.0  # 8 + 10 = 18, so min 22
    assert "center_x" in box
    assert "center_y" in box
    assert "center_z" in box


# ─────────────────────────────────────────────────────────────────────────────
# Manifest-driven path resolution (Codex Milestone 1)
# ─────────────────────────────────────────────────────────────────────────────

# Manifest columns: asset_id  seq_id  gene  stage  asset_type  path  selection_reason  sha256
_MANIFEST_CONTENT = """\
asset_id\tseq_id\tgene\tstage\tasset_type\tpath\tselection_reason\tsha256
MCP_AF2_PDB\tNA23_RS01195\tmcpB\t4\tpdb\t{mcp_pdb}\tTier1 canonical structure\tplaceholder
MCP_P2RANK\tNA23_RS01195\tmcpB\t5\tp2rank_csv\t{mcp_csv}\tTier1 canonical pocket\tplaceholder
"""


def _write_manifest(tmp_path: Path, mcp_pdb: str, mcp_csv: str) -> Path:
    manifest = tmp_path / "aw1_ref_manifest.tsv"
    manifest.write_text(
        _MANIFEST_CONTENT.format(mcp_pdb=mcp_pdb, mcp_csv=mcp_csv),
        encoding="utf-8",
    )
    return manifest


def test_load_aw1_asset_resolves_from_manifest(tmp_path):
    """load_aw1_asset() returns correct path when asset_id is in manifest."""
    # Create a dummy PDB file so the path-exists check passes
    dummy_pdb = tmp_path / "dummy.pdb"
    dummy_pdb.write_text("ATOM dummy\n", encoding="utf-8")

    manifest = _write_manifest(tmp_path, mcp_pdb=str(dummy_pdb), mcp_csv="unused.csv")
    resolved = load_aw1_asset("MCP_AF2_PDB", manifest)
    assert resolved == dummy_pdb.resolve()


def test_load_aw1_asset_raises_for_unknown_id(tmp_path):
    """load_aw1_asset() raises KeyError for unknown asset_id."""
    dummy_pdb = tmp_path / "dummy.pdb"
    dummy_pdb.write_text("ATOM dummy\n", encoding="utf-8")
    manifest = _write_manifest(tmp_path, mcp_pdb=str(dummy_pdb), mcp_csv="unused.csv")

    with pytest.raises(KeyError, match="NONEXISTENT"):
        load_aw1_asset("NONEXISTENT", manifest)


def test_load_aw1_asset_raises_for_missing_manifest(tmp_path):
    """load_aw1_asset() raises FileNotFoundError when manifest file absent."""
    with pytest.raises(FileNotFoundError):
        load_aw1_asset("MCP_AF2_PDB", tmp_path / "no_such_manifest.tsv")


def test_load_aw1_asset_raises_when_path_missing_from_disk(tmp_path):
    """load_aw1_asset() raises FileNotFoundError when the asset file itself is gone."""
    manifest = _write_manifest(
        tmp_path,
        mcp_pdb=str(tmp_path / "ghost.pdb"),   # does NOT exist
        mcp_csv="unused.csv",
    )
    with pytest.raises(FileNotFoundError):
        load_aw1_asset("MCP_AF2_PDB", manifest)


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_load_existing_pdb_manifest_takes_priority(tmp_path):
    """load_existing_pdb() uses manifest path when both manifest and direct path given."""
    manifest = _write_manifest(
        tmp_path,
        mcp_pdb=str(MCP_PDB),
        mcp_csv=str(MCP_PRED_CSV) if PRED_CSV_AVAILABLE else "unused.csv",
    )
    # Pass a nonexistent fallback path — manifest should win
    sp = load_existing_pdb(
        pdb_path="/nonexistent/fallback.pdb",
        seq_id="MCP",
        manifest_path=manifest,
        asset_id="MCP_AF2_PDB",
    )
    assert sp.seq_id == "MCP"
    assert sp.mean_plddt > 0.0
    assert "structure_MCP.pdb" in sp.pdb_path


@pytest.mark.skipif(not ASSETS_AVAILABLE, reason="AW1_ref assets not found")
def test_load_existing_pdb_falls_back_when_manifest_absent(tmp_path):
    """load_existing_pdb() falls back to direct path if manifest not found."""
    sp = load_existing_pdb(
        pdb_path=str(MCP_PDB),
        seq_id="MCP",
        manifest_path=tmp_path / "no_manifest.tsv",  # does not exist
        asset_id="MCP_AF2_PDB",
    )
    assert sp.mean_plddt > 0.0


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_load_existing_predictions_csv_manifest_takes_priority(tmp_path):
    """load_existing_predictions_csv() resolves path via manifest when provided."""
    manifest = _write_manifest(
        tmp_path,
        mcp_pdb=str(MCP_PDB) if ASSETS_AVAILABLE else "unused.pdb",
        mcp_csv=str(MCP_PRED_CSV),
    )
    pockets = load_existing_predictions_csv(
        predictions_csv="/nonexistent/fallback.csv",
        seq_id="MCP",
        manifest_path=manifest,
        asset_id="MCP_P2RANK",
    )
    assert len(pockets) > 0
    assert pockets[0].seq_id == "MCP"


@pytest.mark.skipif(not PRED_CSV_AVAILABLE, reason="AW1_ref predictions CSV not found")
def test_load_existing_predictions_csv_falls_back_when_manifest_absent(tmp_path):
    """load_existing_predictions_csv() falls back to direct path if manifest missing."""
    pockets = load_existing_predictions_csv(
        predictions_csv=str(MCP_PRED_CSV),
        seq_id="MCP",
        manifest_path=tmp_path / "no_manifest.tsv",
        asset_id="MCP_P2RANK",
    )
    assert len(pockets) > 0


@pytest.mark.skipif(
    not (ASSETS_AVAILABLE and PRED_CSV_AVAILABLE),
    reason="AW1_ref assets not found",
)
def test_manifest_driven_full_pipeline(tmp_path):
    """End-to-end: manifest -> structure + pockets -> docking box (Codex Milestone 1)."""
    manifest = _write_manifest(
        tmp_path,
        mcp_pdb=str(MCP_PDB),
        mcp_csv=str(MCP_PRED_CSV),
    )

    sp = load_existing_pdb(
        pdb_path="/nonexistent/fallback.pdb",
        seq_id="MCP",
        manifest_path=manifest,
        asset_id="MCP_AF2_PDB",
    )
    assert sp.mean_plddt > 0.0

    pockets = load_existing_predictions_csv(
        predictions_csv="/nonexistent/fallback.csv",
        seq_id="MCP",
        manifest_path=manifest,
        asset_id="MCP_P2RANK",
    )
    assert len(pockets) > 0

    box = pocket_to_box(pockets[0], ligand_max_dim=8.0)
    assert box["size_x"] >= 22.0
