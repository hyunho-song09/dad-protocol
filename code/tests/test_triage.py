"""
Unit tests for dad.core.triage — Phase B3 validation.

Tests use mock topology results (no real DeepTMHMM/Phobius calls) so they
run without network access or external tool installation.  Mock topologies
are derived from:
  - Phobius literature-based expectations for MCP (nTM=2), Crp (nTM=0), RbsB (SP+nTM=0)
  - Ground truth from AW1_ref/MCP/NA23_RS0_others/NA23_RS06805/phobius_results.png
    (NA23_RS06805: nTM=6, no SP → FLAG with loop extraction or EXCLUDE)

Run with:
    cd d:\\project\\experiment\\DAD\\06_Report\\Mr_Pipeline
    py -3 -m pytest ../Mr_Bio/tests/test_triage.py -v

Results table is printed at end of each test as a summary row.
"""

import sys
import textwrap
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Add dad package to path for import when run from Mr_Bio/tests/
sys.path.insert(0, str(Path(__file__).parents[2] / "Mr_Pipeline"))

from dad.core.triage import (
    FunctionalClass,
    PriorityTier,
    TriageConfig,
    TriageRecord,
    TriageStatus,
    TopologyResult,
    apply_triage_rules,
    assign_priority,
    clip_signal_peptide,
    extract_dock_eligible_region,
    status_to_verdict,
    write_triage_report,
    write_pass_fasta,
)

# ---------------------------------------------------------------------------
# Real sequences from 01_Tier1_input/proteins/
# (concat multi-line FASTA, strip whitespace)
# ---------------------------------------------------------------------------

SEQ_MCP = (
    "MRNMSIFMKVMVIVLILALGMIVIGVYSTFALRNNITEKTMQNLKALAENSGENLVSFIEQ"
    "HTKLIDMLSRDANVMGVYKNEHEEDVWMKKLFNTVLKSYPDVMYVYVGLKDKRMYLIPETE"
    "LPEGYDPTIRPWYQAAVAKPGQVIITEPYADASTGQLVVTVAKTIQTDEGIVGVVALDFDI"
    "SKLSEKLMTKGKELGYLNAVVSKEGNIIMHSDKTLVGKNVANEEFFKKWMSGDESGVFGYT"
    "LNGVKRISGYKRLSNGWIFATLVLEKSLMKDVNRATFINISITIIAILLGIVVALFVTRGF"
    "VVKPISFLVEAAERIANGDLTTKIDYTAKDEIGKLATALSKMVNSLREIALNIERDSATVK"
    "QEASQVAAVSEEVSATIEELTAQVDNVNTNVNNASAAIEEMTSGIEEVAASAQNVAHASQK"
    "LSEEAQKVSNLANEGQKAILSISDVIVQTRQKADATFRIVEQLSESAKNIGEIVDTINSIA"
    "EQTNLLALNAAIEAARAGEAGRGFAVVADEIRKLAEESKQATQNIANILRGIVDSSMKASE"
    "ATKETVEIVNKAYSESDLVKSQFEQILQSIVRMSQMTENLAASAQEQSAAAEEMSSAMDSA"
    "SKSMVSVVEQMNEVTMAIKQQADAISNVARTVENLDNIAEKLVETVRRFKI"
)

SEQ_CRP = (
    "METFVLPAGTKIYDYMEVPRYLYYLIDGIVETKLFNRTMRIETGALGEWALFNLPSQEQVV"
    "SVSEVTLFAIKPDEIYNFEHTGKALKNIIPSVAARLLLIDSELTESQEIPTYVGPDRMRFF"
    "NKVHPGAYKITDSIFQDMLQVKSLYAAGYYKEAYDVIVKLMPQTINEELRKEIMIWHTLLS"
    "MILEPEKADVHFRRLSPKDYSEHLSYLYLTSFYKGGEKQEILEIYMKAGLHLPSETIVTLE"
    "GEVATEAFLVLKGYLKAVKLFEDREVLMSIVGPGEFVGEGAIMNSKTRMATLYSISPTDII"
    "PLSTESIEKAALTNPGFILKICESQLRRIMQVKQLLEIKSQGNQIQRTIMAINYFKPIFGK"
    "AKVSVRDIAYLVDVNVERVIEEVKRMGYKIGIDGSIGV"
)

SEQ_RBSB = (
    "MKRLVLVLAVVLLAVYAFSFKVGLSLSTLNNPFFVTLRDGAMAAAKELGITLLVVDAQDKP"
    "AKQLNDIEDLIQKKVDLIIINPTDSAAIVPAIEAANKAKIPVITVDRAAAGGQVVVHIASD"
    "NVAGGAMAAKFIAEQLKGKGKVVMLVGIPGTSAARDRGTGFKTELKKYPGIQLVAEQVANF"
    "NRAEGMRVMENILQAYPDIDAVFAQNDEMALGAIEAIRAAKKLGKIIVVGFDAIPDALEAV"
    "KKGEMAATVAQQPYLMGQLAVQKAYEYLKTKTIFIPVELELVKK"
)

# NA23_RS06805 — 6-TM MCP (ground truth from phobius_results.png)
# Phobius output: TRANSMEM at 12-35, 47-65, 77-102, 114-136, 157-174, 194-217
# Total length from FASTA: not in our 3 test proteins but included for ground truth test
SEQ_RS06805_FASTA = Path(
    "D:/project/experiment/DAD/AW1_ref/MCP/NA23_RS0_others/NA23_RS06805/"
    "NA23_RS06805 methyl-accepting chemo.fasta"
)

# ---------------------------------------------------------------------------
# Mock topology helpers
# ---------------------------------------------------------------------------

def _mock_mcp_topology() -> TopologyResult:
    """
    Mock DeepTMHMM-equivalent topology for NA23_RS01195 (MCP, nTM=2).

    Biology: MCP has N-terminal periplasmic sensory domain (~160 aa),
    then TM1 (~161-183), short cytoplasmic linker (~184-209),
    TM2 (~210-229), then long cytoplasmic signalling domain.
    The N-terminal periplasmic domain is the dock-eligible region.
    """
    return TopologyResult(
        orf_id="NA23_RS01195",
        n_tm=2,
        has_signal_peptide=False,
        sp_cleavage_site=None,
        tm_regions=[(161, 183), (210, 229)],
        topology_string="O" * 160 + "M" * 23 + "i" * 26 + "M" * 20 + "i" * 432,
        tool_used="mock",
    )


def _mock_crp_topology() -> TopologyResult:
    """
    Mock topology for NA23_RS08105 (Crp/Fnr, nTM=0).
    Crp/Fnr is a soluble cytoplasmic transcriptional regulator.
    No signal peptide, no TM helices.
    """
    return TopologyResult(
        orf_id="NA23_RS08105",
        n_tm=0,
        has_signal_peptide=False,
        sp_cleavage_site=None,
        tm_regions=[],
        topology_string="i" * len(SEQ_CRP),
        tool_used="mock",
    )


def _mock_rbsb_topology() -> TopologyResult:
    """
    Mock topology for NA23_RS00870 (RbsB SBP, SP+nTM=0).
    RbsB is exported to the periplasm via a signal peptide.
    No TM helices in mature form.
    SP cleavage site ~22 aa (MKRLVLVLAVVLLAVYAFSF|KVGL...).
    """
    return TopologyResult(
        orf_id="NA23_RS00870",
        n_tm=0,
        has_signal_peptide=True,
        sp_cleavage_site=22,   # 1-based: SP occupies residues 1-22
        tm_regions=[],
        topology_string="n" * 22 + "o" * (len(SEQ_RBSB) - 22),
        tool_used="mock",
    )


def _mock_rs06805_topology(seq_len: int) -> TopologyResult:
    """
    Mock topology for NA23_RS06805 based on Phobius ground truth.
    phobius_results.png shows: 6 TM helices, no SP.
    TM spans (from Phobius): 12-35, 47-65, 77-102, 114-136, 157-174, 194-217
    Large cytoplasmic C-terminal domain from 218 to end (~354 aa).
    """
    return TopologyResult(
        orf_id="NA23_RS06805",
        n_tm=6,
        has_signal_peptide=False,
        sp_cleavage_site=None,
        tm_regions=[(12, 35), (47, 65), (77, 102), (114, 136), (157, 174), (194, 217)],
        topology_string="mock_6TM",
        tool_used="mock",
    )


# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------

def _run_rule(
    orf_id: str,
    seq: str,
    topo: TopologyResult,
    func_class: FunctionalClass,
    config: Optional[TriageConfig] = None,
) -> TriageRecord:
    if config is None:
        config = TriageConfig()
    return apply_triage_rules(orf_id, seq, topo, func_class, config)


def _print_result_row(label: str, record: TriageRecord) -> None:
    print(
        f"\n{'='*70}\n"
        f"ORF     : {label}\n"
        f"status  : {record.status.value}\n"
        f"nTM     : {record.n_tm}\n"
        f"SP      : {record.has_signal_peptide}\n"
        f"dock    : [{record.dock_region_start}:{record.dock_region_end}]"
        f"  ({record.dock_length} aa)\n"
        f"priority: {record.priority_tier.value}\n"
        f"func    : {record.functional_class.value}\n"
        f"notes   : {record.notes}"
    )


# ---------------------------------------------------------------------------
# Tests — R1: length filter
# ---------------------------------------------------------------------------

class TestR1LengthFilter:

    def test_short_sequence_excluded(self):
        seq = "ACDEFGHIKLMNPQRSTVWY" * 2  # 40 aa < 50
        topo = TopologyResult(
            orf_id="short", n_tm=0, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[],
            topology_string="i" * 40, tool_used="mock",
        )
        record = _run_rule("short", seq, topo, FunctionalClass.HYPOTHETICAL)
        _print_result_row("SHORT (40 aa)", record)
        assert record.status == TriageStatus.EXCLUDE
        assert record.dock_length == 0
        assert "R1" in record.notes

    def test_boundary_49aa_excluded(self):
        seq = "A" * 49
        topo = TopologyResult(
            orf_id="x", n_tm=0, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[], topology_string="i"*49,
            tool_used="mock",
        )
        record = _run_rule("x", seq, topo, FunctionalClass.HYPOTHETICAL)
        assert record.status == TriageStatus.EXCLUDE

    def test_boundary_50aa_passes(self):
        seq = "A" * 50
        topo = TopologyResult(
            orf_id="x", n_tm=0, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[], topology_string="i"*50,
            tool_used="mock",
        )
        record = _run_rule("x", seq, topo, FunctionalClass.OTHER)
        assert record.status != TriageStatus.EXCLUDE or "R1" not in record.notes


# ---------------------------------------------------------------------------
# Tests — R2: signal peptide clipping
# ---------------------------------------------------------------------------

class TestR2SignalPeptide:

    def test_clip_signal_peptide_basic(self):
        seq = "MKFLLL" * 4 + "AREALPROTEIN" * 10  # 24 + 120 = 144 aa
        mature, start, end = clip_signal_peptide(seq, cleavage_site=24)
        assert start == 25              # 1-based
        assert end == len(seq)
        assert len(mature) == len(seq) - 24
        assert mature[0] == seq[24]    # 0-based index 24 = residue 25

    def test_rbsb_sp_clip(self):
        """RbsB: SP ~22 aa, mature form should be ~266 aa."""
        topo = _mock_rbsb_topology()
        record = _run_rule(
            "NA23_RS00870", SEQ_RBSB, topo, FunctionalClass.SBP
        )
        _print_result_row("RbsB (SBP, SP+nTM=0)", record)
        assert record.status == TriageStatus.PASS_CLIPPED
        assert record.has_signal_peptide is True
        assert record.dock_region_start == 23  # 1-based residue after SP
        assert record.dock_length == len(SEQ_RBSB) - 22
        assert record.priority_tier == PriorityTier.HIGH
        assert "R2" in record.notes


# ---------------------------------------------------------------------------
# Tests — R3 + R4: TM topology and dock-region extraction
# ---------------------------------------------------------------------------

class TestR3TmTopology:

    def test_crp_soluble_pass(self):
        """Crp/Fnr: nTM=0, no SP → full sequence PASS, HIGH priority."""
        topo = _mock_crp_topology()
        record = _run_rule(
            "NA23_RS08105", SEQ_CRP, topo, FunctionalClass.REGULATOR_CRP
        )
        _print_result_row("Crp (CRP/FNR, nTM=0)", record)
        assert record.status == TriageStatus.PASS
        assert record.n_tm == 0
        assert record.has_signal_peptide is False
        assert record.dock_region_start == 1
        assert record.dock_region_end == len(SEQ_CRP)
        assert record.dock_length == len(SEQ_CRP)
        assert record.priority_tier == PriorityTier.HIGH
        assert "R3" in record.notes
        assert "nTM=0" in record.notes

    def test_mcp_two_tm_domain_extraction(self):
        """
        MCP: nTM=2 → extract periplasmic N-terminal domain [1:160] = 160 aa.
        Periplasmic domain >= 60 aa → PASS_CLIPPED, HIGH priority.
        """
        topo = _mock_mcp_topology()
        record = _run_rule(
            "NA23_RS01195", SEQ_MCP, topo, FunctionalClass.SENSOR_MCP
        )
        _print_result_row("MCP (nTM=2, periplasmic domain)", record)
        assert record.status == TriageStatus.PASS_CLIPPED
        assert record.n_tm == 2
        assert record.dock_region_start == 1
        assert record.dock_region_end == 160
        assert record.dock_length == 160
        assert record.priority_tier == PriorityTier.HIGH
        assert "R3" in record.notes

    def test_polytopic_6tm_flag_or_exclude(self):
        """
        NA23_RS06805 (Phobius ground truth): nTM=6.
        C-terminal cytoplasmic domain (218-571) = 354 aa >= 60 aa → FLAG.
        """
        # Read real FASTA if available; otherwise use synthetic sequence
        if SEQ_RS06805_FASTA.exists():
            lines = SEQ_RS06805_FASTA.read_text().splitlines()
            seq = "".join(l.strip() for l in lines if not l.startswith(">"))
        else:
            seq = "A" * 571  # synthetic; length from Phobius output (max position 571)

        topo = _mock_rs06805_topology(len(seq))
        record = _run_rule(
            "NA23_RS06805", seq, topo, FunctionalClass.SENSOR_MCP
        )
        _print_result_row("NA23_RS06805 (Phobius ground truth, nTM=6)", record)
        # C-terminal domain 218-571 = 354 aa >= min_loop=60 → extracted, FLAG or PASS_CLIPPED
        assert record.status in (TriageStatus.FLAG, TriageStatus.PASS_CLIPPED, TriageStatus.PASS)
        assert record.n_tm == 6
        assert record.dock_length > 0   # domain should be extracted

    def test_polytopic_7tm_exclude(self):
        """nTM >= 7 → always EXCLUDE."""
        seq = "A" * 400
        tm_regions = [(i*50 + 1, i*50 + 20) for i in range(7)]
        topo = TopologyResult(
            orf_id="poly7", n_tm=7, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=tm_regions,
            topology_string="mock_7TM", tool_used="mock",
        )
        record = _run_rule("poly7", seq, topo, FunctionalClass.OTHER)
        _print_result_row("Polytopic 7TM → EXCLUDE", record)
        assert record.status == TriageStatus.EXCLUDE
        assert "R3" in record.notes
        assert "polytopic" in record.notes.lower()

    def test_single_tm_no_extracellular_domain_exclude(self):
        """nTM=1 with only a 20 aa non-TM region → EXCLUDE (< 60 aa loop)."""
        seq = "A" * 50   # total 50 aa
        # TM at 21-40, leaving only 20 aa N-terminal and 10 aa C-terminal
        topo = TopologyResult(
            orf_id="short_ecto", n_tm=1, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[(21, 40)],
            topology_string="O"*20 + "M"*20 + "i"*10, tool_used="mock",
        )
        record = _run_rule("short_ecto", seq, topo, FunctionalClass.OTHER)
        _print_result_row("Single-TM, short ecto → EXCLUDE", record)
        assert record.status == TriageStatus.EXCLUDE

    def test_multipass_3tm_loop_extracted(self):
        """nTM=3-6: extracts longest non-TM loop >= 60 aa."""
        seq = "A" * 300
        # TM at 10-30, 50-70, 100-120 → non-TM loops: [1-9]=9, [31-49]=19, [71-99]=29, [121-300]=180
        topo = TopologyResult(
            orf_id="multi3", n_tm=3, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[(10, 30), (50, 70), (100, 120)],
            topology_string="mock_3TM", tool_used="mock",
        )
        record = _run_rule("multi3", seq, topo, FunctionalClass.ENZYME)
        _print_result_row("Multi-pass 3TM → loop extraction", record)
        # [121-300] = 180 aa → largest loop, extracted
        assert record.status != TriageStatus.EXCLUDE or "R3" not in record.notes
        if record.status not in (TriageStatus.EXCLUDE,):
            assert record.dock_length == 180


# ---------------------------------------------------------------------------
# Tests — R4: dock-region minimum length
# ---------------------------------------------------------------------------

class TestR4DockLength:

    def test_dock_region_exactly_50_passes(self):
        """Dock region exactly 50 aa at boundary should pass R4."""
        seq = "A" * 200
        # TM at 51-170 → N-terminal 50 aa (barely passes R3 min_loop=60? No, 50 < 60)
        # Use nTM=0 + SP to create a 50 aa mature form
        topo = TopologyResult(
            orf_id="x", n_tm=0, has_signal_peptide=True, sp_cleavage_site=150,
            tm_regions=[], topology_string="n"*150 + "o"*50, tool_used="mock",
        )
        record = _run_rule("x", seq, topo, FunctionalClass.OTHER)
        assert record.dock_length == 50
        assert record.status != TriageStatus.EXCLUDE or "R4" not in record.notes

    def test_dock_region_49_excluded(self):
        """Dock region 49 aa → R4 EXCLUDE."""
        seq = "A" * 200
        topo = TopologyResult(
            orf_id="x", n_tm=0, has_signal_peptide=True, sp_cleavage_site=151,
            tm_regions=[], topology_string="n"*151 + "o"*49, tool_used="mock",
        )
        record = _run_rule("x", seq, topo, FunctionalClass.OTHER)
        _print_result_row("Dock region 49 aa → R4 EXCLUDE", record)
        assert record.status == TriageStatus.EXCLUDE
        assert "R4" in record.notes


# ---------------------------------------------------------------------------
# Tests — R5: functional class priority
# ---------------------------------------------------------------------------

class TestR5FunctionalClass:

    def test_sensor_mcp_high_priority(self):
        record = apply_triage_rules(
            "test_mcp", SEQ_MCP,
            _mock_mcp_topology(), FunctionalClass.SENSOR_MCP, TriageConfig()
        )
        assert record.priority_tier == PriorityTier.HIGH

    def test_regulator_crp_high_priority(self):
        record = apply_triage_rules(
            "test_crp", SEQ_CRP,
            _mock_crp_topology(), FunctionalClass.REGULATOR_CRP, TriageConfig()
        )
        assert record.priority_tier == PriorityTier.HIGH

    def test_sbp_high_priority(self):
        record = apply_triage_rules(
            "test_rbsb", SEQ_RBSB,
            _mock_rbsb_topology(), FunctionalClass.SBP, TriageConfig()
        )
        assert record.priority_tier == PriorityTier.HIGH

    def test_hypothetical_flag(self):
        """Hypothetical protein → FLAG, LOW priority."""
        seq = "A" * 100
        topo = TopologyResult(
            orf_id="hypo", n_tm=0, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[], topology_string="i"*100,
            tool_used="mock",
        )
        record = _run_rule("hypo", seq, topo, FunctionalClass.HYPOTHETICAL)
        _print_result_row("Hypothetical → FLAG LOW", record)
        assert record.status == TriageStatus.FLAG
        assert record.priority_tier == PriorityTier.LOW

    def test_assign_priority_mapping(self):
        assert assign_priority(FunctionalClass.SENSOR_MCP) == PriorityTier.HIGH
        assert assign_priority(FunctionalClass.SENSOR_HK) == PriorityTier.HIGH
        assert assign_priority(FunctionalClass.REGULATOR_CRP) == PriorityTier.HIGH
        assert assign_priority(FunctionalClass.REGULATOR_RR) == PriorityTier.HIGH
        assert assign_priority(FunctionalClass.SBP) == PriorityTier.HIGH
        assert assign_priority(FunctionalClass.ENZYME) == PriorityTier.NORMAL
        assert assign_priority(FunctionalClass.OTHER) == PriorityTier.NORMAL
        assert assign_priority(FunctionalClass.HYPOTHETICAL) == PriorityTier.LOW
        assert assign_priority(FunctionalClass.TRANSPORTER_POLYTOPIC) == PriorityTier.EXCLUDED


# ---------------------------------------------------------------------------
# Tests — extract_dock_eligible_region
# ---------------------------------------------------------------------------

class TestExtractDockRegion:

    def test_mcp_nterm_preferred(self):
        """For nTM=2 MCP, N-terminal periplasmic domain is preferred."""
        seq = "A" * 661
        topo = _mock_mcp_topology()
        config = TriageConfig()
        result = extract_dock_eligible_region(seq, topo, config)
        assert result is not None
        domain_seq, start, end = result
        assert start == 1
        assert end == 160
        assert len(domain_seq) == 160

    def test_no_tm_returns_full(self):
        seq = "A" * 200
        topo = TopologyResult(
            orf_id="x", n_tm=0, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[], topology_string="i"*200,
            tool_used="mock",
        )
        result = extract_dock_eligible_region(seq, topo, TriageConfig())
        assert result is not None
        _, start, end = result
        assert start == 1 and end == 200

    def test_no_eligible_region_returns_none(self):
        seq = "A" * 100
        # TM regions cover all but short gaps (each < 60 aa)
        topo = TopologyResult(
            orf_id="x", n_tm=2, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[(1, 50), (55, 100)],
            topology_string="mock", tool_used="mock",
        )
        result = extract_dock_eligible_region(seq, topo, TriageConfig())
        assert result is None  # only 4 aa between TM1 and TM2


# ---------------------------------------------------------------------------
# Integration test — all 3 Tier1 proteins together
# ---------------------------------------------------------------------------

class TestTier1IntegrationTable:

    def test_three_proteins_expected_outcomes(self):
        """
        Integration test: run all 3 Tier 1 proteins and verify expected
        status, dock region, and priority.  Prints a results table.
        """
        test_cases = [
            ("NA23_RS01195", SEQ_MCP, _mock_mcp_topology(), FunctionalClass.SENSOR_MCP,
             TriageStatus.PASS_CLIPPED, 1, 160, PriorityTier.HIGH),
            ("NA23_RS08105", SEQ_CRP, _mock_crp_topology(), FunctionalClass.REGULATOR_CRP,
             TriageStatus.PASS, 1, len(SEQ_CRP), PriorityTier.HIGH),
            ("NA23_RS00870", SEQ_RBSB, _mock_rbsb_topology(), FunctionalClass.SBP,
             TriageStatus.PASS_CLIPPED, 23, len(SEQ_RBSB), PriorityTier.HIGH),
        ]

        print("\n\n" + "="*80)
        print("TIER 1 TRIAGE INTEGRATION RESULTS")
        print("="*80)
        print(
            f"{'ORF':<18} {'Total':>6} {'Status':<14} "
            f"{'Dock start':>10} {'Dock end':>8} {'Dock len':>8} {'Priority':<8} {'Notes (truncated)'}"
        )
        print("-"*80)

        config = TriageConfig()
        for (orf_id, seq, topo, func, exp_status, exp_start, exp_end, exp_priority) in test_cases:
            record = apply_triage_rules(orf_id, seq, topo, func, config)
            print(
                f"{orf_id:<18} {len(seq):>6} {record.status.value:<14} "
                f"{record.dock_region_start:>10} {record.dock_region_end:>8} "
                f"{record.dock_length:>8} {record.priority_tier.value:<8} "
                f"{record.notes[:40]}..."
            )
            assert record.status == exp_status, (
                f"{orf_id}: expected {exp_status}, got {record.status}"
            )
            assert record.dock_region_start == exp_start, (
                f"{orf_id}: expected dock_start={exp_start}, got {record.dock_region_start}"
            )
            assert record.dock_region_end == exp_end, (
                f"{orf_id}: expected dock_end={exp_end}, got {record.dock_region_end}"
            )
            assert record.priority_tier == exp_priority, (
                f"{orf_id}: expected priority={exp_priority}, got {record.priority_tier}"
            )

        print("="*80)
        print("All 3 Tier 1 proteins passed integration assertions.\n")


# ---------------------------------------------------------------------------
# Tests — report writers
# ---------------------------------------------------------------------------

class TestReportWriters:

    def test_write_triage_report_tsv(self, tmp_path):
        config = TriageConfig()
        records = [
            apply_triage_rules(
                "NA23_RS01195", SEQ_MCP, _mock_mcp_topology(),
                FunctionalClass.SENSOR_MCP, config
            ),
            apply_triage_rules(
                "NA23_RS08105", SEQ_CRP, _mock_crp_topology(),
                FunctionalClass.REGULATOR_CRP, config
            ),
            apply_triage_rules(
                "NA23_RS00870", SEQ_RBSB, _mock_rbsb_topology(),
                FunctionalClass.SBP, config
            ),
        ]
        tsv_path = tmp_path / "triage_report.tsv"
        out = write_triage_report(records, tsv_path)
        assert out.exists()
        lines = out.read_text(encoding="utf-8").splitlines()
        assert lines[0].startswith("orf_id\t")          # header
        assert len(lines) == 4                            # 1 header + 3 data rows
        assert "NA23_RS01195" in lines[1]
        assert "NA23_RS08105" in lines[2]
        assert "NA23_RS00870" in lines[3]

        # Verify both verdict and triage_status columns present in header
        header_cols = lines[0].split("\t")
        assert "verdict" in header_cols
        assert "triage_status" in header_cols
        verdict_idx = header_cols.index("verdict")
        triage_status_idx = header_cols.index("triage_status")

        # MCP: PASS_CLIPPED -> verdict=accept, triage_status=PASS_CLIPPED
        row1 = lines[1].split("\t")
        assert row1[verdict_idx] == "accept"
        assert row1[triage_status_idx] == "PASS_CLIPPED"

        # Crp: PASS -> verdict=accept, triage_status=PASS
        row2 = lines[2].split("\t")
        assert row2[verdict_idx] == "accept"
        assert row2[triage_status_idx] == "PASS"

        print(f"\nTSV output ({tsv_path}):\n" + "\n".join(lines[:5]))

    def test_write_pass_fasta(self, tmp_path):
        config = TriageConfig()
        records = [
            apply_triage_rules(
                "NA23_RS01195", SEQ_MCP, _mock_mcp_topology(),
                FunctionalClass.SENSOR_MCP, config
            ),
            apply_triage_rules(
                "NA23_RS08105", SEQ_CRP, _mock_crp_topology(),
                FunctionalClass.REGULATOR_CRP, config
            ),
        ]
        fasta_path = tmp_path / "triage_pass.fasta"
        out = write_pass_fasta(records, fasta_path)
        assert out.exists()
        text = out.read_text()
        assert ">NA23_RS01195" in text
        assert ">NA23_RS08105" in text
        # Verify FASTA header format
        header_line = [l for l in text.splitlines() if l.startswith(">NA23_RS01195")][0]
        parts = header_line.lstrip(">").split("|")
        assert parts[1] == "PASS_CLIPPED"
        assert parts[3] == "HIGH"
        print(f"\nPass FASTA ({fasta_path}):\n{text[:400]}...")


# ---------------------------------------------------------------------------
# Verdict mapping tests (Codex feedback: Milestone 2 schema fix)
# ---------------------------------------------------------------------------

class TestVerdictMapping:
    """
    Verify status_to_verdict() maps all four TriageStatus values to the
    correct dad.io.TriageDecision verdict strings:
        PASS         -> "accept"
        PASS_CLIPPED -> "accept"
        FLAG         -> "downrank"
        EXCLUDE      -> "exclude"
    """

    def test_pass_maps_to_accept(self):
        assert status_to_verdict(TriageStatus.PASS) == "accept"

    def test_pass_clipped_maps_to_accept(self):
        assert status_to_verdict(TriageStatus.PASS_CLIPPED) == "accept"

    def test_flag_maps_to_downrank(self):
        assert status_to_verdict(TriageStatus.FLAG) == "downrank"

    def test_exclude_maps_to_exclude(self):
        assert status_to_verdict(TriageStatus.EXCLUDE) == "exclude"

    def test_all_statuses_covered(self):
        """Every TriageStatus value must have a mapping (no KeyError)."""
        for status in TriageStatus:
            verdict = status_to_verdict(status)
            assert verdict in {"accept", "downrank", "exclude"}, (
                f"Unexpected verdict {verdict!r} for status {status}"
            )

    def test_tsv_verdict_column_for_flag(self, tmp_path):
        """FLAG status in TSV must have verdict=downrank, triage_status=FLAG."""
        seq = "A" * 100
        topo = TopologyResult(
            orf_id="hypo", n_tm=0, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[], topology_string="i"*100,
            tool_used="mock",
        )
        record = apply_triage_rules("hypo", seq, topo, FunctionalClass.HYPOTHETICAL, TriageConfig())
        assert record.status == TriageStatus.FLAG

        tsv_path = tmp_path / "flag_test.tsv"
        write_triage_report([record], tsv_path)
        lines = tsv_path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        row = lines[1].split("\t")
        assert row[header.index("verdict")] == "downrank"
        assert row[header.index("triage_status")] == "FLAG"

    def test_tsv_verdict_column_for_exclude(self, tmp_path):
        """EXCLUDE status in TSV must have verdict=exclude, triage_status=EXCLUDE."""
        seq = "A" * 20   # < 50 aa -> R1 EXCLUDE
        topo = TopologyResult(
            orf_id="tiny", n_tm=0, has_signal_peptide=False,
            sp_cleavage_site=None, tm_regions=[], topology_string="i"*20,
            tool_used="mock",
        )
        record = apply_triage_rules("tiny", seq, topo, FunctionalClass.OTHER, TriageConfig())
        assert record.status == TriageStatus.EXCLUDE

        tsv_path = tmp_path / "exclude_test.tsv"
        write_triage_report([record], tsv_path)
        lines = tsv_path.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        row = lines[1].split("\t")
        assert row[header.index("verdict")] == "exclude"
        assert row[header.index("triage_status")] == "EXCLUDE"


# ---------------------------------------------------------------------------
# Config and edge-case tests
# ---------------------------------------------------------------------------

class TestConfig:

    def test_invalid_triage_tool_raises(self):
        from dad.core.triage import get_topology
        import pytest
        config = TriageConfig(triage_tool="UnknownTool")
        with pytest.raises(ValueError, match="Unknown triage_tool"):
            get_topology({"x": "ACDEF" * 20}, config)

    def test_default_config_values(self):
        config = TriageConfig()
        assert config.triage_tool == "DeepTMHMM"
        assert config.min_total_length == 50
        assert config.max_tm_exclude == 7
        assert config.min_loop_length_multipass == 60
        assert config.min_dock_region_length == 50


# ---------------------------------------------------------------------------
# Main — print summary table when run directly
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import logging
    logging.basicConfig(level=logging.INFO)

    config = TriageConfig()
    proteins = {
        "NA23_RS01195 (MCP)": (SEQ_MCP, _mock_mcp_topology(), FunctionalClass.SENSOR_MCP),
        "NA23_RS08105 (Crp)": (SEQ_CRP, _mock_crp_topology(), FunctionalClass.REGULATOR_CRP),
        "NA23_RS00870 (RbsB)": (SEQ_RBSB, _mock_rbsb_topology(), FunctionalClass.SBP),
    }

    print("\n" + "="*90)
    print("DAD Stage 3 Triage - Phase B3 Unit Test Summary (mock topology)")
    print("="*90)
    print(
        f"{'Protein':<25} {'Len':>5} {'nTM':>4} {'SP':>3} "
        f"{'Status':<14} {'Dock':>10} {'Len':>6} {'Priority':<8}"
    )
    print("-"*90)

    for label, (seq, topo, func) in proteins.items():
        r = apply_triage_rules(label.split()[0], seq, topo, func, config)
        dock_range = f"{r.dock_region_start}-{r.dock_region_end}" if r.dock_length else "EXCL"
        print(
            f"{label:<25} {len(seq):>5} {r.n_tm:>4} {str(r.has_signal_peptide):>3} "
            f"{r.status.value:<14} {dock_range:>10} {r.dock_length:>6} "
            f"{r.priority_tier.value:<8}"
        )

    print("="*90)
    print("\nRun with: py -3 -m pytest 06_Report/Mr_Bio/tests/test_triage.py -v")
