"""
dad.core.triage — Stage 3: Pre-docking biological triage.

Primary novelty of the DAD protocol: automatically selects dock-competent ORFs
from a translation pool by predicting membrane topology (DeepTMHMM or Phobius)
and applying five decision rules (R1–R5).

Decision rules
--------------
R1  Total length >= 50 aa
R2  Clip signal peptide if present (dock mature form)
R3  nTM=0 → PASS; nTM=1-2 → extract periplasmic domain >= 60 aa;
    nTM=3-6 → extract longest non-TM loop >= 60 aa; nTM >= 7 → EXCLUDE
R4  Dock-eligible region >= 50 aa
R5  Functional class (HMMER/Pfam) → priority tier (HIGH/NORMAL/LOW/EXCLUDED)

Tool selection
--------------
Default : DeepTMHMM (MIT licence, containerisable via biolib)
Opt-in  : Phobius (requires local install; academic free)
Config  : TriageConfig.triage_tool = "DeepTMHMM" | "Phobius"
          or config.yaml: triage_tool: deeptmhmm|phobius

Implementation owner: Mr_Bio (Phase B).
Biological rationale: 06_Report/Mr_Bio/rationale.md
"""

from __future__ import annotations

import csv
import logging
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

logger = logging.getLogger("dad.triage")

# ---------------------------------------------------------------------------
# Enumerations and data classes
# ---------------------------------------------------------------------------

class TriageStatus(str, Enum):
    PASS = "PASS"
    PASS_CLIPPED = "PASS_CLIPPED"
    FLAG = "FLAG"
    EXCLUDE = "EXCLUDE"


class FunctionalClass(str, Enum):
    SENSOR_MCP = "SENSOR_MCP"
    SENSOR_HK = "SENSOR_HK"
    REGULATOR_CRP = "REGULATOR_CRP"
    REGULATOR_RR = "REGULATOR_RR"
    SBP = "SBP"
    ENZYME = "ENZYME"
    TRANSPORTER_POLYTOPIC = "TRANSPORTER_POLYTOPIC"
    HYPOTHETICAL = "HYPOTHETICAL"
    OTHER = "OTHER"


class PriorityTier(str, Enum):
    HIGH = "HIGH"
    NORMAL = "NORMAL"
    LOW = "LOW"
    EXCLUDED = "EXCLUDED"


@dataclass
class TopologyResult:
    """Raw topology prediction output for one ORF."""
    orf_id: str
    n_tm: int
    has_signal_peptide: bool
    sp_cleavage_site: Optional[int]      # 1-based; None if no SP
    tm_regions: List[Tuple[int, int]]    # 1-based (start, end) spans
    topology_string: str                  # raw per-residue annotation
    tool_used: str                        # "DeepTMHMM" | "Phobius" | "mock"


@dataclass
class TriageRecord:
    """Complete triage decision for one ORF (internal + TSV export)."""
    orf_id: str
    status: TriageStatus
    n_tm: int
    has_signal_peptide: bool
    dock_seq: str
    dock_region_start: int               # 1-based
    dock_region_end: int                 # 1-based
    dock_length: int
    functional_class: FunctionalClass
    priority_tier: PriorityTier
    notes: str
    tool_used: str
    source_tier: str = "Tier1"


# ---------------------------------------------------------------------------
# Configuration (zero-parameter defaults)
# ---------------------------------------------------------------------------

@dataclass
class TriageConfig:
    """
    All tunable thresholds.  Defaults are biologically justified in rationale.md.
    Override via CLI or config.yaml — end-users should not need to change these.
    """
    triage_tool: str = "DeepTMHMM"      # or "Phobius"
    phobius_exe: Optional[str] = None   # path to phobius.pl; None → web API
    phobius_url: str = "https://phobius.sbc.su.se/cgi-bin/predict.pl"
    deeptmhmm_exe: Optional[str] = None # path to deeptmhmm binary; None → biolib

    min_total_length: int = 50
    max_tm_exclude: int = 7
    min_loop_length_multipass: int = 60
    min_dock_region_length: int = 50

    # Pfam accessions → priority (R5 rule)
    high_priority_pfam: Tuple[str, ...] = (
        "PF00015",  # MCP_signal
        "PF00325",  # Chemotax_reg
        "PF00027",  # cNMP_binding (CRP/FNR)
        "PF00072",  # Response_reg
        "PF01547",  # SBP_bac_1
        "PF00496",  # SBP_bac_3 (RbsB)
    )
    exclude_pfam: Tuple[str, ...] = (
        "PF00005",  # ABC_tran ATPase — cytoplasmic, no metabolite pocket
    )

    # Pfam accession → FunctionalClass mapping
    pfam_to_class: Dict[str, str] = field(default_factory=lambda: {
        "PF00015": "SENSOR_MCP",
        "PF00325": "SENSOR_HK",
        "PF00027": "REGULATOR_CRP",
        "PF00072": "REGULATOR_RR",
        "PF01547": "SBP",
        "PF00496": "SBP",
        "PF00005": "TRANSPORTER_POLYTOPIC",
        "PF00069": "ENZYME",
        "PF00106": "ENZYME",
    })


# ---------------------------------------------------------------------------
# Topology prediction — DeepTMHMM
# ---------------------------------------------------------------------------

def _write_fasta(sequences: Dict[str, str], path: Path) -> None:
    with open(path, "w") as fh:
        for orf_id, seq in sequences.items():
            fh.write(f">{orf_id}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def _parse_deeptmhmm_gff3(
    gff3_text: str,
    sequences: Dict[str, str],
) -> Dict[str, TopologyResult]:
    """
    Parse DeepTMHMM GFF3 output (produced by `deeptmhmm --gff3` or biolib).

    DeepTMHMM GFF3 feature types:
        TMhelix   — transmembrane helix
        signal    — signal peptide
        inside    — cytoplasmic segment
        outside   — periplasmic / extracellular segment
        Beta      — beta-barrel TM strand

    Example lines:
        NA23_RS01195  DeepTMHMM  signal    1   22  .  .  .  Note=SP
        NA23_RS01195  DeepTMHMM  outside   23  160 .  .  .  .
        NA23_RS01195  DeepTMHMM  TMhelix   161 183 .  .  .  .
    """
    results: Dict[str, TopologyResult] = {}
    current: Dict[str, dict] = {}

    for line in gff3_text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        seq_id, _src, feat, start, end = parts[0], parts[1], parts[2], int(parts[3]), int(parts[4])
        if seq_id not in current:
            current[seq_id] = {
                "n_tm": 0,
                "has_sp": False,
                "sp_cleavage": None,
                "tm_regions": [],
                "topology_chars": [],
            }
        entry = current[seq_id]
        feat_lower = feat.lower()
        if feat_lower in ("tmhelix", "beta"):
            entry["n_tm"] += 1
            entry["tm_regions"].append((start, end))
            entry["topology_chars"].append(f"M{start}-{end}")
        elif feat_lower == "signal":
            entry["has_sp"] = True
            entry["sp_cleavage"] = end  # cleavage after last SP residue
            entry["topology_chars"].append(f"SP{start}-{end}")
        elif feat_lower == "outside":
            entry["topology_chars"].append(f"O{start}-{end}")
        elif feat_lower == "inside":
            entry["topology_chars"].append(f"i{start}-{end}")

    for seq_id, entry in current.items():
        seq = sequences.get(seq_id, "")
        results[seq_id] = TopologyResult(
            orf_id=seq_id,
            n_tm=entry["n_tm"],
            has_signal_peptide=entry["has_sp"],
            sp_cleavage_site=entry["sp_cleavage"],
            tm_regions=sorted(entry["tm_regions"]),
            topology_string="|".join(entry["topology_chars"]),
            tool_used="DeepTMHMM",
        )
    return results


def run_deeptmhmm(
    sequences: Dict[str, str],
    config: TriageConfig,
) -> Dict[str, TopologyResult]:
    """
    Run DeepTMHMM on {orf_id: sequence}.

    Tries (in order):
        1. Local binary at config.deeptmhmm_exe
        2. biolib Python API (lazy import; requires `pip install biolib`)

    Raises
    ------
    RuntimeError
        If neither the local binary nor biolib is available.
    """
    try:
        import biolib  # type: ignore  # lazy import — optional dependency
    except ImportError:
        biolib = None

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_fasta = Path(tmpdir) / "input.fasta"
        _write_fasta(sequences, tmp_fasta)

        # --- Option 1: local binary ---
        if config.deeptmhmm_exe:
            exe = Path(config.deeptmhmm_exe)
            if not exe.exists():
                raise FileNotFoundError(
                    f"deeptmhmm binary not found at {config.deeptmhmm_exe}. "
                    "Install DeepTMHMM or set config.deeptmhmm_exe = None to use biolib."
                )
            out_gff3 = Path(tmpdir) / "results.gff3"
            result = subprocess.run(
                [str(exe), "--fasta", str(tmp_fasta), "--gff3", str(out_gff3)],
                capture_output=True, text=True, timeout=300,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"DeepTMHMM subprocess failed (exit {result.returncode}):\n"
                    f"{result.stderr}"
                )
            gff3_text = out_gff3.read_text()
            return _parse_deeptmhmm_gff3(gff3_text, sequences)

        # --- Option 2: biolib Python API ---
        if biolib is None:
            raise RuntimeError(
                "DeepTMHMM is not available. Install biolib (`pip install biolib`) "
                "or provide config.deeptmhmm_exe pointing to a local DeepTMHMM binary. "
                "Alternatively, set config.triage_tool = 'Phobius' with a local install."
            )

        logger.info("Running DeepTMHMM via biolib API (requires network)...")
        model = biolib.load("DTU/DeepTMHMM")
        job = model.cli(args=f"--fasta {tmp_fasta}")
        job.save_files(tmpdir)

        # biolib saves results as results.gff3 in the output dir
        gff3_candidates = list(Path(tmpdir).glob("**/*.gff3"))
        if not gff3_candidates:
            raise RuntimeError(
                "DeepTMHMM biolib run completed but no GFF3 output was found. "
                "Check biolib job logs."
            )
        gff3_text = gff3_candidates[0].read_text()
        return _parse_deeptmhmm_gff3(gff3_text, sequences)


# ---------------------------------------------------------------------------
# Topology prediction — Phobius
# ---------------------------------------------------------------------------

def _parse_phobius_short(
    text: str,
    sequences: Dict[str, str],
) -> Dict[str, TopologyResult]:
    """
    Parse Phobius short-format output.

    Short-format line:
        SEQUENCE <id>  <length>  <nTM>  <SP>  <prediction>

    prediction chars:
        i = intracellular  o = extracellular/periplasmic
        H = TM helix       n = signal peptide region
    """
    results: Dict[str, TopologyResult] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("SEQENCE") or line.startswith("//"):
            continue
        parts = line.split()
        if len(parts) < 6 or parts[0] != "SEQUENCE":
            continue
        orf_id = parts[1]
        n_tm_raw = int(parts[3])
        has_sp = parts[4] == "Y"
        topo_str = parts[5] if len(parts) > 5 else ""

        # Parse SP cleavage site: 'n' region ends at cleavage
        sp_cleavage: Optional[int] = None
        if has_sp:
            n_count = 0
            for ch in topo_str:
                if ch == "n":
                    n_count += 1
                else:
                    break
            sp_cleavage = n_count if n_count > 0 else None

        # Parse TM regions from topology string
        tm_regions: List[Tuple[int, int]] = []
        in_tm = False
        tm_start = 0
        for idx, ch in enumerate(topo_str, start=1):
            if ch == "H" and not in_tm:
                in_tm = True
                tm_start = idx
            elif ch != "H" and in_tm:
                in_tm = False
                tm_regions.append((tm_start, idx - 1))
        if in_tm:
            tm_regions.append((tm_start, len(topo_str)))

        results[orf_id] = TopologyResult(
            orf_id=orf_id,
            n_tm=n_tm_raw,
            has_signal_peptide=has_sp,
            sp_cleavage_site=sp_cleavage,
            tm_regions=tm_regions,
            topology_string=topo_str,
            tool_used="Phobius",
        )
    return results


def run_phobius(
    sequences: Dict[str, str],
    config: TriageConfig,
) -> Dict[str, TopologyResult]:
    """
    Run Phobius on {orf_id: sequence}.

    Tries (in order):
        1. Local executable at config.phobius_exe (phobius.pl)
        2. Phobius web API (POST to config.phobius_url; rate-limited to 1 req/s)

    Raises
    ------
    RuntimeError  If neither option succeeds.
    FileNotFoundError  If phobius_exe is set but not found.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_fasta = Path(tmpdir) / "input.fasta"
        _write_fasta(sequences, tmp_fasta)

        # --- Option 1: local Phobius perl script ---
        if config.phobius_exe:
            exe = Path(config.phobius_exe)
            if not exe.exists():
                raise FileNotFoundError(
                    f"Phobius executable not found at {config.phobius_exe}. "
                    "Install Phobius from https://phobius.sbc.su.se/ (academic licence)."
                )
            result = subprocess.run(
                ["perl", str(exe), "-short", str(tmp_fasta)],
                capture_output=True, text=True, timeout=120,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"Phobius subprocess failed (exit {result.returncode}):\n"
                    f"{result.stderr}"
                )
            return _parse_phobius_short(result.stdout, sequences)

        # --- Option 2: Phobius web API ---
        try:
            import requests  # type: ignore  # lazy import
        except ImportError:
            raise RuntimeError(
                "Phobius web API requires the 'requests' package (`pip install requests`). "
                "Alternatively, set config.phobius_exe to a local Phobius installation."
            )

        logger.info("Running Phobius via web API (rate-limited)...")
        all_output_lines: List[str] = []
        for orf_id, seq in sequences.items():
            fasta_str = f">{orf_id}\n{seq}\n"
            resp = requests.post(
                config.phobius_url,
                data={"format": "short", "sequence": fasta_str},
                timeout=30,
            )
            if resp.status_code != 200:
                logger.warning(
                    f"Phobius web API returned HTTP {resp.status_code} for {orf_id}; skipping."
                )
                continue
            all_output_lines.extend(resp.text.splitlines())
            time.sleep(1.0)  # respect rate limit

        return _parse_phobius_short("\n".join(all_output_lines), sequences)


# ---------------------------------------------------------------------------
# Topology dispatcher
# ---------------------------------------------------------------------------

def get_topology(
    sequences: Dict[str, str],
    config: TriageConfig,
) -> Dict[str, TopologyResult]:
    """Dispatch topology prediction to the configured tool."""
    tool = config.triage_tool.lower()
    if tool == "deeptmhmm":
        return run_deeptmhmm(sequences, config)
    elif tool == "phobius":
        return run_phobius(sequences, config)
    else:
        raise ValueError(
            f"Unknown triage_tool: {config.triage_tool!r}. "
            "Valid options: 'DeepTMHMM', 'Phobius'."
        )


# ---------------------------------------------------------------------------
# Signal peptide clipping (R2)
# ---------------------------------------------------------------------------

def clip_signal_peptide(seq: str, cleavage_site: int) -> Tuple[str, int, int]:
    """
    Remove signal peptide.  cleavage_site is 1-based (SP ends here).
    Returns (mature_seq, new_start_1based, new_end_1based).
    """
    mature_seq = seq[cleavage_site:]   # 0-based Python slice
    new_start = cleavage_site + 1      # 1-based
    new_end = len(seq)
    return mature_seq, new_start, new_end


# ---------------------------------------------------------------------------
# Non-TM domain extraction (R3)
# ---------------------------------------------------------------------------

def extract_dock_eligible_region(
    seq: str,
    topo: TopologyResult,
    config: TriageConfig,
) -> Optional[Tuple[str, int, int]]:
    """
    Extract the longest non-TM segment that meets the minimum loop length.

    For nTM=2 (MCP-like): prefers the N-terminal region before TM1
    (typically the periplasmic sensory domain).

    Returns (domain_seq, start_1based, end_1based) or None.
    """
    if not topo.tm_regions:
        return seq, 1, len(seq)

    sorted_tm = sorted(topo.tm_regions)

    # Build non-TM intervals
    non_tm: List[Tuple[int, int]] = []
    prev_end = 0
    for tm_start, tm_end in sorted_tm:
        if tm_start - 1 > prev_end:
            non_tm.append((prev_end + 1, tm_start - 1))
        prev_end = tm_end
    if prev_end < len(seq):
        non_tm.append((prev_end + 1, len(seq)))

    # For MCP-like nTM=2: prefer N-terminal region (before TM1)
    if topo.n_tm == 2 and non_tm:
        n_term_start, n_term_end = non_tm[0]
        n_term_len = n_term_end - n_term_start + 1
        if n_term_len >= config.min_loop_length_multipass:
            domain_seq = seq[n_term_start - 1:n_term_end]
            return domain_seq, n_term_start, n_term_end

    # General case: longest non-TM interval >= threshold
    best: Optional[Tuple[int, int]] = None
    best_len = 0
    for start, end in non_tm:
        length = end - start + 1
        if length >= config.min_loop_length_multipass and length > best_len:
            best = (start, end)
            best_len = length

    if best is None:
        return None

    start, end = best
    domain_seq = seq[start - 1:end]
    return domain_seq, start, end


# ---------------------------------------------------------------------------
# Functional annotation (R5) — HMMER/Pfam
# ---------------------------------------------------------------------------

def annotate_functional_class(
    orf_id: str,
    seq: str,
    config: TriageConfig,
) -> FunctionalClass:
    """
    Run hmmscan against Pfam-A and map top hit to FunctionalClass.

    Falls back to HYPOTHETICAL if HMMER is not installed or Pfam-A is
    not available.  This is intentionally permissive — downstream stages
    can always override.

    HMMER must be in PATH; Pfam-A.hmm must be accessible.
    """
    try:
        # Check hmmscan availability (lazy)
        result = subprocess.run(
            ["hmmscan", "--version"],
            capture_output=True, text=True, timeout=5,
        )
        if result.returncode != 0:
            raise FileNotFoundError
    except (FileNotFoundError, subprocess.TimeoutExpired):
        logger.debug(
            f"{orf_id}: hmmscan not found; functional class defaults to HYPOTHETICAL. "
            "Install HMMER + Pfam-A.hmm for R5 annotation."
        )
        return FunctionalClass.HYPOTHETICAL

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_fasta = Path(tmpdir) / "query.fasta"
        tmp_fasta.write_text(f">{orf_id}\n{seq}\n")
        domtbl = Path(tmpdir) / "hits.domtbl"

        # Use --cut_ga (gathering threshold) for Pfam-curated hits only
        cmd = [
            "hmmscan", "--domtblout", str(domtbl),
            "--cut_ga", "--cpu", "1",
            "Pfam-A.hmm", str(tmp_fasta),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            logger.warning(
                f"{orf_id}: hmmscan failed (exit {result.returncode}). "
                "Pfam-A.hmm may not be in the working directory. "
                "Defaulting to HYPOTHETICAL."
            )
            return FunctionalClass.HYPOTHETICAL

        # Parse domtblout: first non-comment, non-empty line with a Pfam accession
        top_pfam: Optional[str] = None
        for line in domtbl.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            cols = line.split()
            if len(cols) >= 4:
                top_pfam = cols[1].split(".")[0]  # strip version suffix
                break

        if top_pfam is None:
            return FunctionalClass.HYPOTHETICAL

        class_name = config.pfam_to_class.get(top_pfam)
        if class_name is None:
            return FunctionalClass.HYPOTHETICAL
        try:
            return FunctionalClass(class_name)
        except ValueError:
            return FunctionalClass.HYPOTHETICAL


def assign_priority(func_class: FunctionalClass) -> PriorityTier:
    _HIGH = {
        FunctionalClass.SENSOR_MCP,
        FunctionalClass.SENSOR_HK,
        FunctionalClass.REGULATOR_CRP,
        FunctionalClass.REGULATOR_RR,
        FunctionalClass.SBP,
    }
    _EXCL = {FunctionalClass.TRANSPORTER_POLYTOPIC}
    _LOW = {FunctionalClass.HYPOTHETICAL}

    if func_class in _EXCL:
        return PriorityTier.EXCLUDED
    if func_class in _HIGH:
        return PriorityTier.HIGH
    if func_class in _LOW:
        return PriorityTier.LOW
    return PriorityTier.NORMAL


# ---------------------------------------------------------------------------
# Core decision tree (R1–R5)
# ---------------------------------------------------------------------------

def apply_triage_rules(
    orf_id: str,
    seq: str,
    topo: TopologyResult,
    func_class: FunctionalClass,
    config: TriageConfig,
    source_tier: str = "Tier1",
) -> TriageRecord:
    """Apply decision rules R1–R5 and return a TriageRecord."""
    notes: List[str] = []
    dock_seq = seq
    dock_start = 1
    dock_end = len(seq)

    # R1 — total length
    if len(seq) < config.min_total_length:
        logger.debug(f"{orf_id}: R1 EXCLUDE — length {len(seq)} < {config.min_total_length}")
        return TriageRecord(
            orf_id=orf_id, status=TriageStatus.EXCLUDE,
            n_tm=topo.n_tm, has_signal_peptide=topo.has_signal_peptide,
            dock_seq="", dock_region_start=0, dock_region_end=0, dock_length=0,
            functional_class=func_class, priority_tier=PriorityTier.EXCLUDED,
            notes=f"R1: length {len(seq)} < {config.min_total_length} aa",
            tool_used=topo.tool_used, source_tier=source_tier,
        )

    # R2 — signal peptide clipping
    if topo.has_signal_peptide and topo.sp_cleavage_site is not None:
        dock_seq, dock_start, dock_end = clip_signal_peptide(seq, topo.sp_cleavage_site)
        notes.append(f"R2: SP clipped at site {topo.sp_cleavage_site}; mature form docked")
        logger.debug(f"{orf_id}: R2 SP clipped → [{dock_start}:{dock_end}]")

    # R3 — TM topology
    n_tm = topo.n_tm

    if n_tm == 0:
        notes.append("R3: nTM=0 — soluble/cytoplasmic; full mature sequence passed")
        logger.debug(f"{orf_id}: R3 PASS soluble (nTM=0)")

    elif 1 <= n_tm <= 2:
        region = extract_dock_eligible_region(dock_seq, topo, config)
        if region is None:
            return TriageRecord(
                orf_id=orf_id, status=TriageStatus.EXCLUDE,
                n_tm=n_tm, has_signal_peptide=topo.has_signal_peptide,
                dock_seq="", dock_region_start=0, dock_region_end=0, dock_length=0,
                functional_class=func_class, priority_tier=PriorityTier.EXCLUDED,
                notes=(f"R3: nTM={n_tm} — no non-TM loop "
                       f">= {config.min_loop_length_multipass} aa; EXCLUDE"),
                tool_used=topo.tool_used, source_tier=source_tier,
            )
        extracted_seq, rel_start, rel_end = region
        dock_seq = extracted_seq
        # Convert relative coords (within potentially-clipped dock_seq) to original space
        abs_start = dock_start + rel_start - 1
        abs_end = abs_start + len(extracted_seq) - 1
        dock_start, dock_end = abs_start, abs_end
        notes.append(
            f"R3: nTM={n_tm} — periplasmic/extracellular domain extracted "
            f"[{dock_start}:{dock_end}]"
        )
        logger.debug(f"{orf_id}: R3 domain extracted [{dock_start}:{dock_end}]")

    elif 3 <= n_tm <= 6:
        region = extract_dock_eligible_region(dock_seq, topo, config)
        if region is None:
            return TriageRecord(
                orf_id=orf_id, status=TriageStatus.FLAG,
                n_tm=n_tm, has_signal_peptide=topo.has_signal_peptide,
                dock_seq=dock_seq, dock_region_start=dock_start, dock_region_end=dock_end,
                dock_length=len(dock_seq),
                functional_class=func_class, priority_tier=PriorityTier.LOW,
                notes=(f"R3: nTM={n_tm} — no loop >= {config.min_loop_length_multipass} aa; "
                       "FLAG for manual review"),
                tool_used=topo.tool_used, source_tier=source_tier,
            )
        extracted_seq, rel_start, rel_end = region
        dock_seq = extracted_seq
        abs_start = dock_start + rel_start - 1
        abs_end = abs_start + len(extracted_seq) - 1
        dock_start, dock_end = abs_start, abs_end
        notes.append(
            f"R3: nTM={n_tm} — largest non-TM loop extracted [{dock_start}:{dock_end}]"
        )
        logger.debug(f"{orf_id}: R3 multi-pass loop extracted [{dock_start}:{dock_end}]")

    else:  # n_tm >= config.max_tm_exclude
        return TriageRecord(
            orf_id=orf_id, status=TriageStatus.EXCLUDE,
            n_tm=n_tm, has_signal_peptide=topo.has_signal_peptide,
            dock_seq="", dock_region_start=0, dock_region_end=0, dock_length=0,
            functional_class=func_class, priority_tier=PriorityTier.EXCLUDED,
            notes=f"R3: nTM={n_tm} >= {config.max_tm_exclude} — polytopic integral membrane; EXCLUDE",
            tool_used=topo.tool_used, source_tier=source_tier,
        )

    # R4 — dock-region length
    dock_length = len(dock_seq)
    if dock_length < config.min_dock_region_length:
        return TriageRecord(
            orf_id=orf_id, status=TriageStatus.EXCLUDE,
            n_tm=n_tm, has_signal_peptide=topo.has_signal_peptide,
            dock_seq="", dock_region_start=0, dock_region_end=0, dock_length=0,
            functional_class=func_class, priority_tier=PriorityTier.EXCLUDED,
            notes=(f"R4: dock-eligible region {dock_length} aa "
                   f"< {config.min_dock_region_length} aa; EXCLUDE"),
            tool_used=topo.tool_used, source_tier=source_tier,
        )

    # R5 — functional class
    priority = assign_priority(func_class)
    if priority == PriorityTier.EXCLUDED:
        return TriageRecord(
            orf_id=orf_id, status=TriageStatus.EXCLUDE,
            n_tm=n_tm, has_signal_peptide=topo.has_signal_peptide,
            dock_seq="", dock_region_start=0, dock_region_end=0, dock_length=0,
            functional_class=func_class, priority_tier=PriorityTier.EXCLUDED,
            notes="R5: TRANSPORTER_POLYTOPIC — cytoplasmic ATPase domain; EXCLUDE",
            tool_used=topo.tool_used, source_tier=source_tier,
        )

    if func_class == FunctionalClass.HYPOTHETICAL:
        status = TriageStatus.FLAG
        notes.append("R5: hypothetical/DUF — FLAG; proceed with low confidence")
    elif topo.has_signal_peptide or n_tm > 0:
        status = TriageStatus.PASS_CLIPPED
        notes.append(f"R5: {func_class.value}; priority {priority.value}")
    else:
        status = TriageStatus.PASS
        notes.append(f"R5: {func_class.value}; priority {priority.value}")

    notes_str = "; ".join(notes)
    logger.info(
        f"{orf_id}: {status.value} | nTM={n_tm} | SP={topo.has_signal_peptide} | "
        f"dock=[{dock_start}:{dock_end}] ({dock_length} aa) | "
        f"{func_class.value} | {priority.value}"
    )
    return TriageRecord(
        orf_id=orf_id, status=status,
        n_tm=n_tm, has_signal_peptide=topo.has_signal_peptide,
        dock_seq=dock_seq,
        dock_region_start=dock_start, dock_region_end=dock_end,
        dock_length=dock_length,
        functional_class=func_class, priority_tier=priority,
        notes=notes_str, tool_used=topo.tool_used, source_tier=source_tier,
    )


# ---------------------------------------------------------------------------
# Report writers
# ---------------------------------------------------------------------------

# Mapping from internal TriageStatus to external I/O verdict (dad.io.TriageDecision)
# PASS/PASS_CLIPPED -> accept, FLAG -> downrank, EXCLUDE -> exclude
_VERDICT_MAP: Dict[TriageStatus, str] = {
    TriageStatus.PASS: "accept",
    TriageStatus.PASS_CLIPPED: "accept",
    TriageStatus.FLAG: "downrank",
    TriageStatus.EXCLUDE: "exclude",
}


def status_to_verdict(status: TriageStatus) -> str:
    """Map internal TriageStatus to dad.io.TriageDecision verdict string."""
    return _VERDICT_MAP[status]


_TSV_COLUMNS = [
    "orf_id", "verdict", "triage_status", "n_tm", "has_signal_peptide",
    "dock_region_start", "dock_region_end", "dock_length",
    "functional_class", "priority_tier", "notes",
    "tool_used", "source_tier",
]


def write_triage_report(
    records: List[TriageRecord],
    output_path: Optional[Path] = None,
) -> Path:
    """Write triage_report.tsv (dock_seq excluded; goes to triage_pass.fasta).

    Two status columns are written for pipeline compatibility:
      verdict       — external I/O schema (accept/downrank/exclude), used by
                      Snakemake helper get_accepted_seq_ids() and dad.io.TriageDecision
      triage_status — internal biological status (PASS/PASS_CLIPPED/FLAG/EXCLUDE),
                      preserved for manuscript and debugging
    """
    if output_path is None:
        output_path = Path("triage_report.tsv")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(_TSV_COLUMNS)
        for r in records:
            writer.writerow([
                r.orf_id, status_to_verdict(r.status), r.status.value,
                r.n_tm, r.has_signal_peptide,
                r.dock_region_start, r.dock_region_end, r.dock_length,
                r.functional_class.value, r.priority_tier.value,
                r.notes, r.tool_used, r.source_tier,
            ])
    logger.info(f"Triage report written to {output_path}")
    return output_path


def write_pass_fasta(
    records: List[TriageRecord],
    output_path: Optional[Path] = None,
) -> Path:
    """
    Write triage_pass.fasta — PASS and PASS_CLIPPED sequences for Stage 4.
    FLAG sequences are included with a warning prefix in the header.

    Header format:
        >{orf_id}|{status}|{dock_start}-{dock_end}|{priority}|{source_tier}
    """
    if output_path is None:
        output_path = Path("triage_pass.fasta")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    _include = {TriageStatus.PASS, TriageStatus.PASS_CLIPPED, TriageStatus.FLAG}
    with open(output_path, "w", encoding="utf-8") as fh:
        for r in records:
            if r.status not in _include or not r.dock_seq:
                continue
            header = (
                f">{r.orf_id}|{r.status.value}|"
                f"{r.dock_region_start}-{r.dock_region_end}|"
                f"{r.priority_tier.value}|{r.source_tier}"
            )
            fh.write(header + "\n")
            for i in range(0, len(r.dock_seq), 60):
                fh.write(r.dock_seq[i:i + 60] + "\n")
    logger.info(f"Pass FASTA written to {output_path}")
    return output_path


# ---------------------------------------------------------------------------
# Main pipeline entry point
# ---------------------------------------------------------------------------

def run_triage(
    sequences: Dict[str, str],
    config: Optional[TriageConfig] = None,
    source_tier: str = "Tier1",
    output_tsv: Optional[Path] = None,
    output_fasta: Optional[Path] = None,
) -> List[TriageRecord]:
    """
    Run full Stage 3 triage pipeline.

    Parameters
    ----------
    sequences   : {orf_id: aa_sequence}
    config      : TriageConfig; uses zero-parameter defaults if None
    source_tier : data provenance label ("Tier1" | "Tier2" | "Tier3")
    output_tsv  : path for triage_report.tsv (default: triage_report.tsv)
    output_fasta: path for triage_pass.fasta (default: triage_pass.fasta)

    Returns
    -------
    List[TriageRecord] — one per ORF including EXCLUDE records.
    """
    if config is None:
        config = TriageConfig()

    logger.info(f"Triage started: {len(sequences)} ORFs via {config.triage_tool}")

    # Step 1: topology prediction (batch)
    topo_results = get_topology(sequences, config)

    records: List[TriageRecord] = []

    # Step 2: per-ORF decision tree
    for orf_id, seq in sequences.items():
        topo = topo_results.get(orf_id)
        if topo is None:
            logger.warning(f"{orf_id}: topology prediction returned no result; treating as FLAG")
            topo = TopologyResult(
                orf_id=orf_id, n_tm=0, has_signal_peptide=False,
                sp_cleavage_site=None, tm_regions=[],
                topology_string="UNKNOWN", tool_used=config.triage_tool,
            )
        func_class = annotate_functional_class(orf_id, seq, config)
        record = apply_triage_rules(orf_id, seq, topo, func_class, config, source_tier)
        records.append(record)

    # Step 3: write outputs
    write_triage_report(records, output_tsv)
    write_pass_fasta(records, output_fasta)

    passed = sum(1 for r in records if r.status in (TriageStatus.PASS, TriageStatus.PASS_CLIPPED))
    flagged = sum(1 for r in records if r.status == TriageStatus.FLAG)
    excluded = sum(1 for r in records if r.status == TriageStatus.EXCLUDE)
    logger.info(
        f"Triage complete: {passed} PASS/CLIPPED | {flagged} FLAG | {excluded} EXCLUDE "
        f"(of {len(records)} total)"
    )
    return records


# ---------------------------------------------------------------------------
# Snakemake rule interface
# ---------------------------------------------------------------------------

def snakemake_triage_main(
    input_fasta: str,
    output_tsv: str,
    output_fasta: str,
    triage_tool: str = "DeepTMHMM",
    source_tier: str = "Tier1",
) -> None:
    """
    Entry point for Snakemake rule 02_triage.smk.

    Snakemake rule usage::

        run:
            from dad.core.triage import snakemake_triage_main
            snakemake_triage_main(
                input.fasta, output.tsv, output.fasta,
                params.triage_tool, params.source_tier,
            )
    """
    try:
        from Bio import SeqIO  # type: ignore
    except ImportError:
        raise RuntimeError(
            "Biopython is required (`pip install biopython`) for FASTA parsing "
            "in the Snakemake interface."
        )

    sequences = {
        str(rec.id): str(rec.seq)
        for rec in SeqIO.parse(input_fasta, "fasta")
    }
    config = TriageConfig(triage_tool=triage_tool)
    run_triage(
        sequences, config, source_tier,
        output_tsv=Path(output_tsv),
        output_fasta=Path(output_fasta),
    )
