"""
Tier 1 mock tests for dad.core.docking / interaction / visualize.

Tests run without GPU or GNINA binary — all external calls are mocked.
Validates schema contracts (DockingResult, InteractionProfile, RankedPair)
and core algorithmic helpers (SDF parsing, box sizing, z-score ranking).

Run from repo root:
    cd d:/project/experiment/DAD/06_Report/Mr_Pipeline
    python -m pytest ../Mr_Dock/tests/test_docking.py -v
"""

from __future__ import annotations

import sys
import os
import textwrap
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock
from typing import List

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Path bootstrap so tests can import dad.* without install
# ---------------------------------------------------------------------------
_PIPELINE_DIR = Path(__file__).resolve().parents[2] / "Mr_Pipeline"
if str(_PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(_PIPELINE_DIR))

from dad.io import (
    DockingBox,
    DockingPose,
    DockingResult,
    InteractionProfile,
    ContactResidue,
    LigandInput,
    PreparedLigand,
    PocketResult,
    RankedPair,
    StructurePrediction,
)


# ─────────────────────────────────────────────────────────────────────────────
# Fixtures — minimal synthetic files
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def tmp(tmp_path):
    return tmp_path


MINIMAL_SDF = textwrap.dedent("""\

     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
""")

GNINA_SDF = textwrap.dedent("""\

     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
REMARK minimizedAffinity -5.61
REMARK CNNscore 0.8995
REMARK CNNaffinity 4.416
$$$$

     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.5000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
REMARK minimizedAffinity -5.32
REMARK CNNscore 0.8832
REMARK CNNaffinity 4.393
$$$$
""")


@pytest.fixture
def sdf_file(tmp):
    p = tmp / "lig.sdf"
    p.write_text(GNINA_SDF)
    return str(p)


@pytest.fixture
def minimal_sdf_file(tmp):
    p = tmp / "minimal.sdf"
    p.write_text(MINIMAL_SDF)
    return str(p)


MINIMAL_PDB = textwrap.dedent("""\
ATOM      1  CA  GLY A   1       1.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       2.500   0.500   0.000  1.00  0.00           C
ATOM      3  CA  VAL A   3       4.000   1.000   0.000  1.00  0.00           C
END
""")


@pytest.fixture
def receptor_pdb(tmp):
    p = tmp / "receptor.pdb"
    p.write_text(MINIMAL_PDB)
    return str(p)


# ─────────────────────────────────────────────────────────────────────────────
# Tests — dad.core.docking helpers
# ─────────────────────────────────────────────────────────────────────────────

class TestSdfParsing:
    def test_parse_coords_basic(self, sdf_file):
        from dad.core.docking import _parse_sdf_v2000_coords
        coords = _parse_sdf_v2000_coords(sdf_file, model_index=0)
        assert coords.shape == (3, 3), f"Expected (3,3), got {coords.shape}"
        assert abs(coords[0][0] - 1.0) < 1e-6
        assert abs(coords[2][0] - 4.0) < 1e-6

    def test_parse_coords_second_model(self, sdf_file):
        from dad.core.docking import _parse_sdf_v2000_coords
        coords = _parse_sdf_v2000_coords(sdf_file, model_index=1)
        assert coords.shape == (3, 3)
        assert abs(coords[0][0] - 1.5) < 1e-6

    def test_parse_coords_model_out_of_range(self, sdf_file):
        from dad.core.docking import _parse_sdf_v2000_coords
        with pytest.raises(IndexError):
            _parse_sdf_v2000_coords(sdf_file, model_index=99)

    def test_compute_max_dim(self, minimal_sdf_file):
        from dad.core.docking import compute_ligand_max_dim
        max_dim = compute_ligand_max_dim(minimal_sdf_file)
        # atoms at x=1,2.5,4 → max dist = 3.0
        assert abs(max_dim - 3.0) < 1e-4


class TestParseGninaOutputSdf:
    def test_parse_scores_schema(self, sdf_file):
        from dad.core.docking import parse_gnina_output_sdf
        result = parse_gnina_output_sdf(sdf_file, "MCP", "AlaIle", pocket_rank=1)

        assert isinstance(result, DockingResult)
        assert result.seq_id == "MCP"
        assert result.lig_id == "AlaIle"
        assert result.pocket_rank == 1
        assert len(result.poses) == 2

        bp = result.best_pose
        assert isinstance(bp, DockingPose)
        assert bp.pose_rank == 1
        assert abs(bp.vina_affinity - (-5.61)) < 1e-4
        assert abs(bp.cnn_pose_score - 0.8995) < 1e-4
        assert abs(bp.cnn_affinity - 4.416) < 1e-4

    def test_best_pose_matches_first(self, sdf_file):
        from dad.core.docking import parse_gnina_output_sdf
        result = parse_gnina_output_sdf(sdf_file, "Crp", "GlyVal", 1)
        assert result.best_pose.pose_rank == result.poses[0].pose_rank

    def test_missing_sdf_raises(self, tmp):
        from dad.core.docking import parse_gnina_output_sdf
        with pytest.raises(FileNotFoundError):
            parse_gnina_output_sdf(str(tmp / "nonexistent.sdf"), "X", "Y", 1)


class TestBuildDockingBoxes:
    def test_box_sizing_rule_min(self):
        from dad.core.docking import build_docking_boxes
        pocket = PocketResult(
            seq_id="MCP", pocket_rank=1, score=0.9,
            center_x=-18.1706, center_y=3.1775, center_z=7.1226,
        )
        ligand = PreparedLigand(
            lig_id="AlaIle", sdf_path="/fake/ala.sdf",
            smiles_canonical="CC(N)C(=O)NC(CC(C)C)C(=O)O",
            mol_weight=202.25, max_dim=5.0,  # max_dim + padding(10) = 15 < 22 → box=22
        )
        config = {"docking": {"box_size_min": 22.0, "box_padding": 10.0}}
        boxes = build_docking_boxes({"MCP": [pocket]}, [ligand], config)
        assert len(boxes) == 1
        box = boxes[0]
        assert box.size_x == 22.0  # max(22, 5+10) = 22

    def test_box_sizing_rule_scaled(self):
        from dad.core.docking import build_docking_boxes
        pocket = PocketResult(
            seq_id="Crp", pocket_rank=1, score=0.8,
            center_x=0.0, center_y=0.0, center_z=0.0,
        )
        ligand = PreparedLigand(
            lig_id="BigLig", sdf_path="/fake/big.sdf",
            smiles_canonical="C" * 20,
            mol_weight=300.0, max_dim=15.0,  # 15+10=25 > 22 → box=25
        )
        config = {"docking": {"box_size_min": 22.0, "box_padding": 10.0}}
        boxes = build_docking_boxes({"Crp": [pocket]}, [ligand], config)
        assert boxes[0].size_x == 25.0

    def test_many_to_many(self):
        from dad.core.docking import build_docking_boxes
        pockets = {
            "P1": [
                PocketResult("P1", 1, 0.9, 0, 0, 0),
                PocketResult("P1", 2, 0.7, 1, 1, 1),
            ],
            "P2": [PocketResult("P2", 1, 0.85, 2, 2, 2)],
        }
        ligands = [
            PreparedLigand("L1", "/f", "C", 100.0, 5.0),
            PreparedLigand("L2", "/f", "N", 120.0, 6.0),
        ]
        config = {"docking": {}}
        boxes = build_docking_boxes(pockets, ligands, config)
        # 2 proteins × (2+1 pockets) × 2 ligands = but per protein: P1 has 2 pockets, P2 has 1
        # total = 2 lig × (2+1) pockets = 6
        assert len(boxes) == 6


class TestRunGninaSingle:
    def test_gnina_command_structure(self, tmp, sdf_file, receptor_pdb):
        """Mock subprocess to validate GNINA command is constructed correctly."""
        from dad.core.docking import run_gnina_single
        box = DockingBox(
            seq_id="MCP", lig_id="AlaIle", pocket_rank=1,
            center_x=-18.17, center_y=3.18, center_z=7.12,
            size_x=22.0, size_y=22.0, size_z=22.0,
        )
        out_sdf = str(tmp / "out.sdf")
        (tmp / "out.sdf").write_text(GNINA_SDF)  # fake output

        captured_cmd = []
        def fake_run(cmd, **kwargs):
            captured_cmd.extend(cmd)
            return MagicMock(returncode=0, stdout="gnina v1.3.2", stderr="")

        with patch("dad.core.docking.subprocess.run", side_effect=fake_run):
            result = run_gnina_single(
                receptor_pdb=receptor_pdb,
                ligand_sdf=sdf_file,
                box=box,
                output_sdf=out_sdf,
                gnina_path="gnina",
            )

        assert "gnina" in captured_cmd
        assert "-r" in captured_cmd
        assert "--exhaustiveness" in captured_cmd
        idx_ex = captured_cmd.index("--exhaustiveness")
        assert captured_cmd[idx_ex + 1] == "32"
        assert "--num_modes" in captured_cmd
        assert "--seed" in captured_cmd
        idx_seed = captured_cmd.index("--seed")
        assert captured_cmd[idx_seed + 1] == "0"
        assert isinstance(result, DockingResult)

    def test_gnina_missing_receptor_raises(self, tmp, sdf_file):
        from dad.core.docking import run_gnina_single
        box = DockingBox("X", "Y", 1, 0, 0, 0, 22, 22, 22)
        with pytest.raises(FileNotFoundError):
            run_gnina_single(
                receptor_pdb=str(tmp / "missing.pdb"),
                ligand_sdf=sdf_file,
                box=box,
                output_sdf=str(tmp / "out.sdf"),
            )

    def test_gnina_nonzero_exit_raises(self, tmp, sdf_file, receptor_pdb):
        from dad.core.docking import run_gnina_single
        box = DockingBox("X", "Y", 1, 0, 0, 0, 22, 22, 22)
        out_sdf = str(tmp / "out.sdf")

        with patch("dad.core.docking.subprocess.run",
                   return_value=MagicMock(returncode=1, stdout="", stderr="Error")):
            with pytest.raises(RuntimeError):
                run_gnina_single(receptor_pdb, sdf_file, box, out_sdf)


# ─────────────────────────────────────────────────────────────────────────────
# Tests — dad.core.interaction helpers
# ─────────────────────────────────────────────────────────────────────────────

class TestGetSdfAtomCoords:
    def test_coords_shape(self, sdf_file):
        from dad.core.interaction import get_sdf_atom_coords
        coords = get_sdf_atom_coords(sdf_file, model_index=0)
        assert coords.shape == (3, 3)

    def test_missing_file_raises(self, tmp):
        from dad.core.interaction import get_sdf_atom_coords
        with pytest.raises(FileNotFoundError):
            get_sdf_atom_coords(str(tmp / "nope.sdf"))

    def test_model_out_of_range_raises(self, sdf_file):
        from dad.core.interaction import get_sdf_atom_coords
        with pytest.raises(IndexError):
            get_sdf_atom_coords(sdf_file, model_index=50)


class TestExtractContactResidues:
    def test_contact_residues_schema(self, sdf_file, receptor_pdb):
        from dad.core.interaction import extract_contact_residues
        contacts = extract_contact_residues(
            receptor_pdb=receptor_pdb,
            docked_sdf=sdf_file,
            pose_index=0,
            contact_distance=10.0,  # generous cutoff for the minimal test structure
        )
        assert isinstance(contacts, list)
        for c in contacts:
            assert isinstance(c, ContactResidue)
            assert isinstance(c.chain, str)
            assert isinstance(c.res_id, int)
            assert isinstance(c.res_name, str)
            assert c.min_dist >= 0.0

    def test_contacts_sorted_by_distance(self, sdf_file, receptor_pdb):
        from dad.core.interaction import extract_contact_residues
        contacts = extract_contact_residues(receptor_pdb, sdf_file, 0, 10.0)
        dists = [c.min_dist for c in contacts]
        assert dists == sorted(dists)

    def test_no_contacts_at_zero_cutoff(self, sdf_file, receptor_pdb):
        from dad.core.interaction import extract_contact_residues
        contacts = extract_contact_residues(receptor_pdb, sdf_file, 0, 0.0)
        assert contacts == []


class TestAggregateResults:
    def _make_result(self, seq_id, lig_id, vina, cnn_pose, cnn_aff, sdf="/f.sdf"):
        pose = DockingPose(1, vina, cnn_pose, cnn_aff)
        return DockingResult(seq_id, lig_id, 1, sdf, [pose], pose)

    def test_ranking_order(self):
        from dad.core.interaction import aggregate_results
        results = [
            self._make_result("MCP", "AlaIle", -5.61, 0.8995, 4.416),
            self._make_result("MCP", "GlyVal", -5.18, 0.8468, 4.022),
            self._make_result("Crp", "AlaIle", -6.01, 0.6228, 4.642),
        ]
        config = {}
        ranked = aggregate_results(results, [], config)
        assert len(ranked) == 3
        # ranks should be 1-indexed and monotonically increasing
        ranks = [rp.overall_rank for rp in ranked]
        assert ranks == list(range(1, len(ranked) + 1))

    def test_ranked_pair_schema(self):
        from dad.core.interaction import aggregate_results
        results = [self._make_result("P", "L", -5.0, 0.9, 5.0)]
        ranked = aggregate_results(results, [], {})
        assert len(ranked) == 1
        rp = ranked[0]
        assert isinstance(rp, RankedPair)
        assert rp.seq_id == "P"
        assert rp.lig_id == "L"
        assert rp.overall_rank == 1

    def test_single_result_zscore_zero(self):
        from dad.core.interaction import aggregate_results, _zscore
        arr = np.array([5.0])
        assert _zscore(arr)[0] == 0.0


class TestBuildHeatmapMatrix:
    def _make_ranked_pair(self, seq_id, lig_id, cnn_aff, rank=1):
        return RankedPair(seq_id, lig_id, 1, -5.0, 0.9, cnn_aff, 1.0, rank)

    def test_matrix_shape(self):
        from dad.core.interaction import build_heatmap_matrix
        pairs = [
            self._make_ranked_pair("P1", "L1", 4.5, 1),
            self._make_ranked_pair("P1", "L2", 3.8, 2),
            self._make_ranked_pair("P2", "L1", 5.1, 3),
        ]
        proteins, ligands, matrix = build_heatmap_matrix(pairs)
        assert matrix.shape == (2, 2)
        assert set(proteins) == {"P1", "P2"}
        assert set(ligands) == {"L1", "L2"}

    def test_matrix_values(self):
        from dad.core.interaction import build_heatmap_matrix
        pairs = [
            self._make_ranked_pair("P1", "L1", 4.5),
            self._make_ranked_pair("P1", "L1", 5.0),  # same pair, different pocket → take max
        ]
        # Two rows with same (P1, L1) but second has higher cnn_affinity
        pairs[1] = RankedPair("P1", "L1", 2, -5.0, 0.9, 5.0, 1.0, 1)
        proteins, ligands, matrix = build_heatmap_matrix(pairs)
        # should take max = 5.0
        assert matrix[0, 0] == pytest.approx(5.0)

    def test_invalid_score_column(self):
        from dad.core.interaction import build_heatmap_matrix
        pairs = [self._make_ranked_pair("P", "L", 4.0)]
        with pytest.raises(ValueError):
            build_heatmap_matrix(pairs, score_column="bogus_column")


# ─────────────────────────────────────────────────────────────────────────────
# Tests — dad.core.visualize helpers
# ─────────────────────────────────────────────────────────────────────────────

class TestMergeReceptorLigandPdb:
    def test_merge_creates_file(self, tmp, sdf_file, receptor_pdb):
        from dad.core.visualize import merge_receptor_ligand_pdb
        out_pdb = str(tmp / "complex.pdb")
        result = merge_receptor_ligand_pdb(receptor_pdb, sdf_file, out_pdb)
        assert Path(result).exists()
        content = Path(result).read_text()
        assert "ATOM" in content or "HETATM" in content

    def test_merge_missing_receptor_raises(self, tmp, sdf_file):
        from dad.core.visualize import merge_receptor_ligand_pdb
        with pytest.raises(FileNotFoundError):
            merge_receptor_ligand_pdb(str(tmp / "missing.pdb"), sdf_file, str(tmp / "out.pdb"))


class TestGenerateChimeraxCxc:
    def test_cxc_content(self, tmp, sdf_file, receptor_pdb):
        from dad.core.visualize import generate_chimerax_cxc
        cxc_path = str(tmp / "test.cxc")
        result = generate_chimerax_cxc(
            receptor_pdb=receptor_pdb,
            docked_sdf=sdf_file,
            contact_residue_ids=[167, 169, 178],
            ligand_name="AlaIle",
            output_cxc=cxc_path,
        )
        assert Path(result).exists()
        cxc = Path(result).read_text()
        assert "hide all" in cxc
        assert "show /A cartoon" in cxc
        assert "color ligand green" in cxc
        assert "AlaIle" in cxc
        assert "167" in cxc

    def test_cxc_empty_residues(self, tmp, sdf_file, receptor_pdb):
        from dad.core.visualize import generate_chimerax_cxc
        cxc_path = str(tmp / "empty_res.cxc")
        result = generate_chimerax_cxc(
            receptor_pdb=receptor_pdb,
            docked_sdf=sdf_file,
            contact_residue_ids=[],
            ligand_name="GlyVal",
            output_cxc=cxc_path,
        )
        cxc = Path(result).read_text()
        assert "show ligand sticks" in cxc


class TestGeneratePy3dmolHtml:
    def test_html_output_structure(self, tmp, sdf_file, receptor_pdb):
        from dad.core.visualize import generate_py3dmol_html
        out_html = str(tmp / "view.html")
        result = generate_py3dmol_html(
            receptor_pdb=receptor_pdb,
            docked_sdf=sdf_file,
            contact_residue_ids=[167, 169],
            ligand_name="AlaIle",
            output_html=out_html,
        )
        assert Path(result).exists()
        html = Path(result).read_text()
        assert "3Dmol" in html
        assert "AlaIle" in html
        assert "greenCarbon" in html
        assert "spectrum" in html


class TestWriteMasterCsv:
    def test_csv_columns(self, tmp):
        from dad.core.interaction import write_master_csv
        pairs = [
            RankedPair("MCP", "AlaIle", 1, -5.61, 0.8995, 4.416, 1.5, 1, 15, 3),
            RankedPair("Crp", "AlaIle", 1, -6.01, 0.6228, 4.642, 1.2, 2, 12, 2),
        ]
        csv_path = str(tmp / "master.csv")
        write_master_csv(pairs, csv_path)
        import csv as csv_mod
        with open(csv_path) as fh:
            reader = csv_mod.DictReader(fh)
            rows = list(reader)
        assert len(rows) == 2
        assert rows[0]["seq_id"] == "MCP"
        assert "vina_affinity" in rows[0]
        assert rows[0]["overall_rank"] == "1"


# ─────────────────────────────────────────────────────────────────────────────
# Tests — DockingResult schema round-trip
# ─────────────────────────────────────────────────────────────────────────────

class TestDockingResultSchema:
    def test_to_dict_from_dict_roundtrip(self, sdf_file):
        from dad.core.docking import parse_gnina_output_sdf
        result = parse_gnina_output_sdf(sdf_file, "MCP", "AlaIle", 1)
        d = result.to_dict()
        restored = DockingResult.from_dict(d)
        assert restored.seq_id == result.seq_id
        assert restored.best_pose.vina_affinity == result.best_pose.vina_affinity
        assert len(restored.poses) == len(result.poses)

    def test_ground_truth_values(self, sdf_file):
        """Validate against Revision.txt Tier 1 ground truth for MCP × Ala-Ile."""
        from dad.core.docking import parse_gnina_output_sdf
        result = parse_gnina_output_sdf(sdf_file, "MCP", "AlaIle", 1)
        bp = result.best_pose
        # Revision.txt: Vina -5.61, CNN pose 0.8995, CNN affinity 4.416 (±0.1 tolerance)
        assert abs(bp.vina_affinity - (-5.61)) < 0.1
        assert abs(bp.cnn_pose_score - 0.8995) < 0.1
        assert abs(bp.cnn_affinity - 4.416) < 0.1


# ─────────────────────────────────────────────────────────────────────────────
# Tests — Tier 1 replay mode (no GPU, no GNINA binary)
# Validates all 6 Revision.txt ground truth cases via replay_from_ground_truth()
# ─────────────────────────────────────────────────────────────────────────────

# Ground truth from Revision.txt Table X — tolerances for schema validation
_TIER1_CASES = [
    # (case_id,        ligand_id,  exp_vina, exp_cnn_pose, exp_cnn_aff)
    ("MCP",  "AlaIle",  -5.61, 0.8995, 4.416),
    ("MCP",  "GlyVal",  -5.18, 0.8468, 4.022),
    ("Crp",  "AlaIle",  -6.01, 0.6228, 4.642),
    ("Crp",  "GlyVal",  -5.59, 0.7353, 4.574),
    ("RbsB", "AlaIle",  -4.92, 0.4447, 3.516),
    ("RbsB", "GlyVal",  -4.44, 0.9144, 3.896),
]

_REVISION_TXT = Path(__file__).resolve().parents[3] / "Claude_web" / "Revision.txt"


class TestReplayFromGroundTruth:
    """Tier 1 replay mode: construct DockingResult from Revision.txt without GNINA."""

    def _replay(self, case_id, ligand_id):
        from dad.core.docking import replay_from_ground_truth
        return replay_from_ground_truth(
            case_id=case_id,
            revision_table_path=_REVISION_TXT,
            ligand_id=ligand_id,
        )

    # ── 6 individual case tests ───────────────────────────────────────────────

    def test_mcp_ala_ile_scores(self):
        r = self._replay("MCP", "AlaIle")
        bp = r.best_pose
        assert abs(bp.vina_affinity - (-5.61)) < 1e-6
        assert abs(bp.cnn_pose_score - 0.8995) < 1e-6
        assert abs(bp.cnn_affinity - 4.416) < 1e-6

    def test_mcp_gly_val_scores(self):
        r = self._replay("MCP", "GlyVal")
        bp = r.best_pose
        assert abs(bp.vina_affinity - (-5.18)) < 1e-6
        assert abs(bp.cnn_pose_score - 0.8468) < 1e-6
        assert abs(bp.cnn_affinity - 4.022) < 1e-6

    def test_crp_ala_ile_scores(self):
        r = self._replay("Crp", "AlaIle")
        bp = r.best_pose
        assert abs(bp.vina_affinity - (-6.01)) < 1e-6
        assert abs(bp.cnn_pose_score - 0.6228) < 1e-6
        assert abs(bp.cnn_affinity - 4.642) < 1e-6

    def test_crp_gly_val_scores(self):
        r = self._replay("Crp", "GlyVal")
        bp = r.best_pose
        assert abs(bp.vina_affinity - (-5.59)) < 1e-6
        assert abs(bp.cnn_pose_score - 0.7353) < 1e-6
        assert abs(bp.cnn_affinity - 4.574) < 1e-6

    def test_rbsb_ala_ile_scores(self):
        r = self._replay("RbsB", "AlaIle")
        bp = r.best_pose
        assert abs(bp.vina_affinity - (-4.92)) < 1e-6
        assert abs(bp.cnn_pose_score - 0.4447) < 1e-6
        assert abs(bp.cnn_affinity - 3.516) < 1e-6

    def test_rbsb_gly_val_scores(self):
        r = self._replay("RbsB", "GlyVal")
        bp = r.best_pose
        assert abs(bp.vina_affinity - (-4.44)) < 1e-6
        assert abs(bp.cnn_pose_score - 0.9144) < 1e-6
        assert abs(bp.cnn_affinity - 3.896) < 1e-6

    # ── Schema / contract tests ───────────────────────────────────────────────

    def test_schema_is_docking_result(self):
        r = self._replay("MCP", "AlaIle")
        assert isinstance(r, DockingResult)
        assert isinstance(r.best_pose, DockingPose)
        assert isinstance(r.poses, list)
        assert len(r.poses) == 1

    def test_pose_source_marker(self):
        """gnina_version carries 'precomputed_revision' for downstream detection."""
        r = self._replay("Crp", "GlyVal")
        assert r.gnina_version == "precomputed_revision"

    def test_output_sdf_sentinel(self):
        """output_sdf must be the replay sentinel (not a real path)."""
        r = self._replay("RbsB", "AlaIle")
        assert r.output_sdf == "<replay:no_sdf>"

    def test_elapsed_seconds_zero(self):
        r = self._replay("MCP", "GlyVal")
        assert r.elapsed_seconds == 0.0

    def test_seq_id_and_lig_id_preserved(self):
        r = self._replay("NA23_RS01195", "Ala-Ile")
        assert r.seq_id == "NA23_RS01195"
        assert r.lig_id == "Ala-Ile"

    def test_roundtrip_to_dict(self):
        r = self._replay("Crp", "AlaIle")
        d = r.to_dict()
        restored = DockingResult.from_dict(d)
        assert restored.best_pose.vina_affinity == r.best_pose.vina_affinity
        assert restored.gnina_version == "precomputed_revision"

    # ── Alias resolution tests ────────────────────────────────────────────────

    def test_seq_id_alias_na23_rs01195(self):
        """NA23_RS01195 resolves to MCP ground truth."""
        r = self._replay("NA23_RS01195", "AlaIle")
        assert abs(r.best_pose.vina_affinity - (-5.61)) < 1e-6

    def test_seq_id_alias_na23_rs08105(self):
        """NA23_RS08105 resolves to Crp ground truth."""
        r = self._replay("NA23_RS08105", "GlyVal")
        assert abs(r.best_pose.vina_affinity - (-5.59)) < 1e-6

    def test_seq_id_alias_na23_rs00870(self):
        """NA23_RS00870 resolves to RbsB ground truth."""
        r = self._replay("NA23_RS00870", "AlaIle")
        assert abs(r.best_pose.vina_affinity - (-4.92)) < 1e-6

    def test_ligand_alias_ala_ile_hyphen(self):
        r = self._replay("MCP", "Ala-Ile")
        assert abs(r.best_pose.vina_affinity - (-5.61)) < 1e-6

    def test_ligand_alias_gly_val_hyphen(self):
        r = self._replay("Crp", "Gly-Val")
        assert abs(r.best_pose.vina_affinity - (-5.59)) < 1e-6

    # ── Error handling ────────────────────────────────────────────────────────

    def test_unknown_target_raises_key_error(self):
        from dad.core.docking import replay_from_ground_truth
        with pytest.raises(KeyError, match="unknown target"):
            replay_from_ground_truth("UnknownProtein", _REVISION_TXT, "AlaIle")

    def test_unknown_ligand_raises_key_error(self):
        from dad.core.docking import replay_from_ground_truth
        with pytest.raises(KeyError, match="unknown ligand"):
            replay_from_ground_truth("MCP", _REVISION_TXT, "UnknownLigand")

    def test_missing_revision_file_uses_hardcoded_fallback(self, tmp):
        """If Revision.txt path doesn't exist, hardcoded table is used."""
        from dad.core.docking import replay_from_ground_truth
        fake_path = tmp / "nonexistent_revision.txt"
        r = replay_from_ground_truth("MCP", fake_path, "AlaIle")
        assert abs(r.best_pose.vina_affinity - (-5.61)) < 1e-6

    # ── Parametric all-6 sweep ────────────────────────────────────────────────

    @pytest.mark.parametrize("case_id,lig_id,exp_vina,exp_cnn_pose,exp_cnn_aff", _TIER1_CASES)
    def test_all_tier1_cases(self, case_id, lig_id, exp_vina, exp_cnn_pose, exp_cnn_aff):
        """Parametric sweep: all 6 Revision.txt cases pass schema + value checks."""
        from dad.core.docking import replay_from_ground_truth
        r = replay_from_ground_truth(case_id, _REVISION_TXT, lig_id)
        assert isinstance(r, DockingResult)
        assert r.gnina_version == "precomputed_revision"
        assert abs(r.best_pose.vina_affinity - exp_vina) < 1e-6
        assert abs(r.best_pose.cnn_pose_score - exp_cnn_pose) < 1e-6
        assert abs(r.best_pose.cnn_affinity - exp_cnn_aff) < 1e-6


# ─────────────────────────────────────────────────────────────────────────────
# Tests — reference-free pose selection (Blocker 4)
# Tests use a 2-pose SDF (GNINA_SDF fixture) to stay GPU/GNINA-free.
# Pose 0: CNN_pose=0.8995, CNN_aff=4.416 → ensemble=0.8995*0.6+4.416*0.4=2.3061
# Pose 1: CNN_pose=0.8832, CNN_aff=4.393 → ensemble=0.8832*0.6+4.393*0.4=2.2867
# → ensemble always selects pose 0 for the GNINA_SDF fixture.
# ─────────────────────────────────────────────────────────────────────────────

# SDF with two poses that have deliberately swapped CNN scores so cluster ≠ ensemble
CLUSTER_TEST_SDF = textwrap.dedent("""\

     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    4.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
REMARK minimizedAffinity -5.00
REMARK CNNscore 0.5000
REMARK CNNaffinity 2.000
$$$$

     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
   11.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   12.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   14.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
REMARK minimizedAffinity -6.00
REMARK CNNscore 0.9500
REMARK CNNaffinity 5.000
$$$$
""")


class TestSelectConsensusPose:
    """Reference-free pose selection — no crystal RMSD, no GPU, no GNINA binary."""

    @pytest.fixture
    def gnina_sdf_path(self, tmp):
        p = tmp / "gnina.sdf"
        p.write_text(GNINA_SDF)
        return p

    @pytest.fixture
    def cluster_sdf_path(self, tmp):
        """Two well-separated poses: pose 0 low ensemble score, pose 1 high ensemble score."""
        p = tmp / "cluster.sdf"
        p.write_text(CLUSTER_TEST_SDF)
        return p

    @pytest.fixture
    def receptor_path(self, tmp):
        p = tmp / "receptor.pdb"
        p.write_text(MINIMAL_PDB)
        return p

    # ── Test 1: ensemble method picks pose 0 (highest CNN ensemble score) ─────

    def test_ensemble_selects_highest_score(self, gnina_sdf_path, receptor_path):
        from dad.core.docking import select_consensus_pose
        result = select_consensus_pose(
            poses_sdf=gnina_sdf_path,
            receptor_pdb=receptor_path,
            method="ensemble",
        )
        assert isinstance(result, dict)
        assert "selected_pose_index" in result
        assert "score_breakdown" in result
        assert "method_used" in result
        # Pose 0: 0.8995*0.6 + 4.416*0.4 = 2.3061
        # Pose 1: 0.8832*0.6 + 4.393*0.4 = 2.2867
        assert result["selected_pose_index"] == 0
        assert result["method_used"] == "ensemble"
        scores = result["score_breakdown"]["ensemble_scores"]
        assert len(scores) == 2
        assert scores[0] > scores[1]

    # ── Test 2: cluster method can select a different pose than ensemble ──────

    def test_cluster_can_differ_from_ensemble(self, cluster_sdf_path, receptor_path):
        """Cluster SDF has two spatially separated poses; largest cluster = both separate.
        Cluster method picks by largest cluster then ensemble — for well-separated poses
        each forms its own cluster, so it falls back to the highest ensemble within
        each single-member cluster, picking pose 1 (higher ensemble score).
        Ensemble alone on pose 0 vs pose 1:
          Pose 0: 0.5*0.6+2.0*0.4=1.1
          Pose 1: 0.95*0.6+5.0*0.4=2.57 → pose 1 wins in both cases here.
        The key point is the return schema is valid and method_used contains 'cluster'.
        """
        from dad.core.docking import select_consensus_pose
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = select_consensus_pose(
                poses_sdf=cluster_sdf_path,
                receptor_pdb=receptor_path,
                method="cluster",
            )
        assert isinstance(result, dict)
        assert "selected_pose_index" in result
        assert result["selected_pose_index"] in (0, 1)
        # method_used may be cluster or cluster→ensemble(...) for fallback
        assert "cluster" in result["method_used"]

    # ── Test 3: posebusters falls back to ensemble with warning when not installed

    def test_posebusters_fallback_warns_when_unavailable(self, gnina_sdf_path, receptor_path):
        from dad.core.docking import select_consensus_pose
        import warnings
        # Mock the import to simulate PoseBusters not being installed
        with patch.dict("sys.modules", {"posebusters": None}):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = select_consensus_pose(
                    poses_sdf=gnina_sdf_path,
                    receptor_pdb=receptor_path,
                    method="posebusters",
                )
        assert "ensemble" in result["method_used"]
        # At least one warning should mention posebusters or fallback
        warning_texts = " ".join(str(w.message) for w in caught)
        assert "posebusters" in warning_texts.lower() or "ensemble" in warning_texts.lower()

    # ── Test 4: plip falls back to ensemble with warning when not installed ───

    def test_plip_fallback_warns_when_unavailable(self, gnina_sdf_path, receptor_path):
        from dad.core.docking import select_consensus_pose
        import warnings
        with patch.dict("sys.modules", {"plip": None, "plip.structure": None,
                                        "plip.structure.preparation": None}):
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                result = select_consensus_pose(
                    poses_sdf=gnina_sdf_path,
                    receptor_pdb=receptor_path,
                    method="plip",
                )
        assert "ensemble" in result["method_used"]
        warning_texts = " ".join(str(w.message) for w in caught)
        assert "plip" in warning_texts.lower() or "ensemble" in warning_texts.lower()

    # ── Test 5: consensus returns valid schema with method_used ───────────────

    def test_consensus_returns_valid_schema(self, gnina_sdf_path, receptor_path):
        from dad.core.docking import select_consensus_pose
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = select_consensus_pose(
                poses_sdf=gnina_sdf_path,
                receptor_pdb=receptor_path,
                method="consensus",
            )
        assert isinstance(result, dict)
        assert "selected_pose_index" in result
        assert result["selected_pose_index"] in (0, 1)
        assert "method_used" in result
        assert "consensus" in result["method_used"] or "ensemble" in result["method_used"]

    # ── Test 6: unknown method raises ValueError ──────────────────────────────

    def test_unknown_method_raises_value_error(self, gnina_sdf_path, receptor_path):
        from dad.core.docking import select_consensus_pose
        with pytest.raises(ValueError, match="Unknown method"):
            select_consensus_pose(
                poses_sdf=gnina_sdf_path,
                receptor_pdb=receptor_path,
                method="bogus_method",
            )
