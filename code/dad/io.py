"""
dad.io — Shared I/O schemas for the DAD pipeline.

All inter-stage data contracts are defined here as dataclasses.
Downstream agents (Mr_Struct, Mr_Dock, Mr_Bio) must import from this module
so that type contracts remain consistent across the pipeline.

Design decisions
----------------
- stdlib ``dataclasses`` only: avoids pydantic dependency at dry-run / scaffold stage.
- All paths are ``str`` (not ``pathlib.Path``) for Snakemake string compatibility.
- Optional fields default to ``None``; mandatory fields have no default.
- Serialization: each dataclass exposes ``to_dict()`` / ``from_dict()`` for JSON round-trips.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1 — Input schemas
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ProteinInput:
    """Single ORF translation record from the input FASTA.

    Parameters
    ----------
    seq_id : str
        FASTA header identifier (no ``>`` prefix, truncated at first space).
    sequence : str
        Amino acid sequence (single-letter codes, no gaps).
    source_file : str
        Absolute path to the originating FASTA file.
    description : str, optional
        Full FASTA header description text.
    organism : str, optional
        Organism name hint provided by the user.
    """

    seq_id: str
    sequence: str
    source_file: str
    description: Optional[str] = None
    organism: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.sequence:
            raise ValueError(f"Empty sequence for {self.seq_id}")
        self.sequence = self.sequence.upper().replace(" ", "")

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "ProteinInput":
        return cls(**d)


@dataclass
class LigandInput:
    """Single ligand record from the input SMILES list.

    Parameters
    ----------
    lig_id : str
        Unique ligand identifier.
    smiles : str
        Canonical SMILES string.
    source_file : str
        Absolute path to the originating SMI/CSV file.
    mol_weight : float, optional
        Molecular weight in Da (computed at prep stage).
    """

    lig_id: str
    smiles: str
    source_file: str
    mol_weight: Optional[float] = None

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "LigandInput":
        return cls(**d)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 3 — Triage schemas
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class TriageDecision:
    """Per-ORF biological triage verdict.

    Parameters
    ----------
    seq_id : str
        Identifier matching ``ProteinInput.seq_id``.
    verdict : str
        External I/O verdict: one of ``"accept"``, ``"downrank"``, ``"exclude"``.
        Mapping from biological status: PASS/PASS_CLIPPED -> accept,
        FLAG -> downrank, EXCLUDE -> exclude.
    confidence : float
        Triage confidence score in [0, 1].
    triage_status : str, optional
        Internal biological status from dad.core.triage:
        ``"PASS"``, ``"PASS_CLIPPED"``, ``"FLAG"``, or ``"EXCLUDE"``.
        Preserved for manuscript and debugging; verdict is the pipeline I/O key.
    topology : str, optional
        Predicted topology string (e.g., ``"SP+2TM"``).
    n_tm : int, optional
        Number of predicted transmembrane helices.
    has_signal_peptide : bool, optional
        Whether a signal peptide was predicted.
    cleavage_site : int, optional
        Signal peptide cleavage position (1-based AA index).
    dock_region_start : int, optional
        First residue index of the dock-competent domain.
    dock_region_end : int, optional
        Last residue index of the dock-competent domain.
    pfam_hits : list of str, optional
        Pfam domain accessions found (e.g., ``["PF00672"]``).
    functional_class : str, optional
        High-level functional class (e.g., ``"MCP"``, ``"Crp_FNR"``, ``"SBP"``).
    priority_flag : bool
        True if functional class matches a priority family.
    rationale : str, optional
        Human-readable explanation for the verdict.
    """

    seq_id: str
    verdict: str                             # "accept" | "downrank" | "exclude"
    confidence: float
    triage_status: Optional[str] = None     # "PASS" | "PASS_CLIPPED" | "FLAG" | "EXCLUDE"
    topology: Optional[str] = None
    n_tm: Optional[int] = None
    has_signal_peptide: Optional[bool] = None
    cleavage_site: Optional[int] = None
    dock_region_start: Optional[int] = None
    dock_region_end: Optional[int] = None
    pfam_hits: Optional[List[str]] = field(default_factory=list)
    functional_class: Optional[str] = None
    priority_flag: bool = False
    rationale: Optional[str] = None

    def __post_init__(self) -> None:
        valid_verdicts = {"accept", "downrank", "exclude"}
        if self.verdict not in valid_verdicts:
            raise ValueError(f"verdict must be one of {valid_verdicts}, got '{self.verdict}'")
        if not (0.0 <= self.confidence <= 1.0):
            raise ValueError(f"confidence must be in [0, 1], got {self.confidence}")
        valid_statuses = {"PASS", "PASS_CLIPPED", "FLAG", "EXCLUDE", None}
        if self.triage_status not in valid_statuses:
            raise ValueError(
                f"triage_status must be one of {valid_statuses}, got '{self.triage_status}'"
            )

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "TriageDecision":
        return cls(**d)

    @classmethod
    def from_tsv_row(cls, row: Dict[str, str]) -> "TriageDecision":
        """Construct from a triage_report.tsv row dict (csv.DictReader output).

        Reads both ``verdict`` and ``triage_status`` columns so that pipeline
        helpers can use the external I/O key (verdict) while preserving the
        internal biological status (triage_status).
        """
        orf_id = row.get("orf_id") or row.get("seq_id", "")
        verdict = row.get("verdict", "")
        triage_status = row.get("triage_status") or row.get("status") or None

        # Infer verdict from triage_status if verdict column missing
        if not verdict and triage_status:
            _map = {
                "PASS": "accept", "PASS_CLIPPED": "accept",
                "FLAG": "downrank", "EXCLUDE": "exclude",
            }
            verdict = _map.get(triage_status, "")

        def _int(v: str) -> Optional[int]:
            try:
                return int(v) if v and v != "0" else None
            except (ValueError, TypeError):
                return None

        def _bool(v: str) -> Optional[bool]:
            if v is None:
                return None
            return v.lower() in ("true", "1", "yes")

        return cls(
            seq_id=orf_id,
            verdict=verdict,
            confidence=1.0,
            triage_status=triage_status or None,
            n_tm=_int(row.get("n_tm", "")),
            has_signal_peptide=_bool(row.get("has_signal_peptide")),
            dock_region_start=_int(row.get("dock_region_start", "")),
            dock_region_end=_int(row.get("dock_region_end", "")),
            functional_class=row.get("functional_class") or None,
            rationale=row.get("notes") or None,
        )


# ─────────────────────────────────────────────────────────────────────────────
# Stage 4 — Structure prediction schemas
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class StructurePrediction:
    """Result of ColabFold / AlphaFold2 structure prediction for one ORF.

    Parameters
    ----------
    seq_id : str
        Identifier matching ``ProteinInput.seq_id``.
    pdb_path : str
        Absolute path to the best-ranked PDB file (rank_1).
    mean_plddt : float
        Mean per-residue pLDDT across all residues.
    ptm_score : float, optional
        Predicted TM-score (for complexes).
    plddt_per_residue : list of float, optional
        Per-residue pLDDT values (length == len(sequence)).
    model_type : str, optional
        ColabFold model type used (e.g., ``"alphafold2_ptm"``).
    num_recycles : int, optional
        Number of recycle iterations performed.
    low_confidence_flag : bool
        True if mean_plddt < config threshold (default 70).
    job_dir : str, optional
        Directory containing all ColabFold output files.
    """

    seq_id: str
    pdb_path: str
    mean_plddt: float
    ptm_score: Optional[float] = None
    plddt_per_residue: Optional[List[float]] = field(default_factory=list)
    model_type: Optional[str] = None
    num_recycles: Optional[int] = None
    low_confidence_flag: bool = False
    job_dir: Optional[str] = None

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "StructurePrediction":
        return cls(**d)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 5 — Pocket detection schemas
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class PocketResult:
    """A single P2Rank-predicted binding pocket.

    Parameters
    ----------
    seq_id : str
        Identifier matching ``ProteinInput.seq_id``.
    pocket_rank : int
        P2Rank rank (1 = top-scoring pocket).
    score : float
        P2Rank pocket probability/score.
    center_x : float
        Pocket centroid x-coordinate (Å).
    center_y : float
        Pocket centroid y-coordinate (Å).
    center_z : float
        Pocket centroid z-coordinate (Å).
    residues : list of int, optional
        Residue indices (1-based) lining the pocket.
    surf_atoms : int, optional
        Number of surface atoms in the pocket.
    """

    seq_id: str
    pocket_rank: int
    score: float
    center_x: float
    center_y: float
    center_z: float
    residues: Optional[List[int]] = field(default_factory=list)
    surf_atoms: Optional[int] = None

    @property
    def center(self) -> Tuple[float, float, float]:
        """Return pocket center as (x, y, z) tuple."""
        return (self.center_x, self.center_y, self.center_z)

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "PocketResult":
        return cls(**d)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 6-7 — Ligand preparation schemas
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class PreparedLigand:
    """3D-prepared ligand ready for docking.

    Parameters
    ----------
    lig_id : str
        Identifier matching ``LigandInput.lig_id``.
    sdf_path : str
        Absolute path to the prepared 3D SDF file.
    smiles_canonical : str
        Canonical SMILES after RDKit standardization.
    mol_weight : float
        Molecular weight in Da.
    max_dim : float
        Maximum atom-pair distance in the 3D conformer (Å); used for box sizing.
    n_rotatable_bonds : int, optional
        Number of rotatable bonds.
    formal_charge : int, optional
        Net formal charge at preparation pH.
    """

    lig_id: str
    sdf_path: str
    smiles_canonical: str
    mol_weight: float
    max_dim: float
    n_rotatable_bonds: Optional[int] = None
    formal_charge: Optional[int] = None

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "PreparedLigand":
        return cls(**d)


@dataclass
class DockingBox:
    """Auto-configured GNINA docking search box.

    Box dimensions follow the DAD rule:
        side = max(22 Å, ligand_max_dim + 10 Å)   [from config.yaml]

    Parameters
    ----------
    seq_id : str
        Target protein identifier.
    lig_id : str
        Ligand identifier.
    pocket_rank : int
        Pocket rank this box is centered on.
    center_x, center_y, center_z : float
        Box center coordinates (from P2Rank output).
    size_x, size_y, size_z : float
        Box dimensions in Å.
    """

    seq_id: str
    lig_id: str
    pocket_rank: int
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "DockingBox":
        return cls(**d)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 8 — Docking result schemas
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class DockingPose:
    """A single GNINA docking pose (one of num_modes outputs).

    Parameters
    ----------
    pose_rank : int
        Pose rank within the docking run (1 = best GNINA-ranked pose).
    vina_affinity : float
        Vina affinity score in kcal/mol (more negative = stronger).
    cnn_pose_score : float
        CNN pose score in [0, 1] (higher = more plausible pose geometry).
    cnn_affinity : float
        CNN-predicted binding affinity.
    rmsd_lb : float, optional
        RMSD to best pose (lower bound cluster), if available.
    rmsd_ub : float, optional
        RMSD to best pose (upper bound cluster), if available.
    """

    pose_rank: int
    vina_affinity: float
    cnn_pose_score: float
    cnn_affinity: float
    rmsd_lb: Optional[float] = None
    rmsd_ub: Optional[float] = None

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "DockingPose":
        return cls(**d)


@dataclass
class DockingResult:
    """Full GNINA docking result for one (protein × ligand × pocket) combination.

    Parameters
    ----------
    seq_id : str
        Target protein identifier.
    lig_id : str
        Ligand identifier.
    pocket_rank : int
        Pocket rank used for this docking run.
    output_sdf : str
        Absolute path to the GNINA output SDF (all poses).
    poses : list of DockingPose
        All returned poses sorted by GNINA rank.
    best_pose : DockingPose
        Convenience alias for poses[0] (best-ranked pose).
    gnina_version : str, optional
        GNINA version string (from ``gnina --version``).
    elapsed_seconds : float, optional
        Wall-clock docking time.
    """

    seq_id: str
    lig_id: str
    pocket_rank: int
    output_sdf: str
    poses: List[DockingPose]
    best_pose: DockingPose
    gnina_version: Optional[str] = None
    elapsed_seconds: Optional[float] = None

    def to_dict(self) -> Dict:
        d = asdict(self)
        return d

    @classmethod
    def from_dict(cls, d: Dict) -> "DockingResult":
        poses = [DockingPose.from_dict(p) for p in d.pop("poses")]
        best_pose = DockingPose.from_dict(d.pop("best_pose"))
        return cls(poses=poses, best_pose=best_pose, **d)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 9 — Interaction profiling schemas
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ContactResidue:
    """A protein residue in contact with the docked ligand.

    Parameters
    ----------
    chain : str
        PDB chain identifier.
    res_id : int
        Residue sequence number (1-based).
    res_name : str
        Three-letter residue name (e.g., ``"GLY"``).
    min_dist : float
        Minimum atom-atom distance to any ligand atom (Å).
    interaction_type : str, optional
        Classified interaction type: ``"hbond"``, ``"hydrophobic"``,
        ``"pi_stack"``, ``"ionic"``, ``"contact"`` (generic).
    """

    chain: str
    res_id: int
    res_name: str
    min_dist: float
    interaction_type: Optional[str] = "contact"

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "ContactResidue":
        return cls(**d)


@dataclass
class InteractionProfile:
    """Full interaction profile for one docked pose.

    Parameters
    ----------
    seq_id : str
        Target protein identifier.
    lig_id : str
        Ligand identifier.
    pocket_rank : int
        Pocket rank.
    pose_rank : int
        Docking pose rank (typically 1).
    contact_residues : list of ContactResidue
        All residues within contact_distance of the ligand.
    n_hbonds : int, optional
        Number of hydrogen bonds (from PLIP).
    n_hydrophobic : int, optional
        Number of hydrophobic contacts (from PLIP).
    plip_xml_path : str, optional
        Path to PLIP XML output for detailed records.
    """

    seq_id: str
    lig_id: str
    pocket_rank: int
    pose_rank: int
    contact_residues: List[ContactResidue]
    n_hbonds: Optional[int] = None
    n_hydrophobic: Optional[int] = None
    plip_xml_path: Optional[str] = None

    def to_dict(self) -> Dict:
        d = asdict(self)
        return d

    @classmethod
    def from_dict(cls, d: Dict) -> "InteractionProfile":
        contacts = [ContactResidue.from_dict(r) for r in d.pop("contact_residues")]
        return cls(contact_residues=contacts, **d)


# ─────────────────────────────────────────────────────────────────────────────
# Stage 10 — Ranking schemas
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class RankedPair:
    """One row in the master ranking table (protein × ligand × pocket).

    Parameters
    ----------
    seq_id : str
        Target protein identifier.
    lig_id : str
        Ligand identifier.
    pocket_rank : int
        Pocket rank used.
    vina_affinity : float
        Best-pose Vina affinity (kcal/mol).
    cnn_pose_score : float
        Best-pose CNN pose score.
    cnn_affinity : float
        Best-pose CNN affinity.
    composite_rank_score : float
        Composite z-score: z(vina)*w1 + z(cnn_pose)*w2 + z(cnn_affinity)*w3.
    overall_rank : int
        Rank across all (protein × ligand × pocket) combinations (1 = best).
    n_contacts : int, optional
        Number of contact residues.
    n_hbonds : int, optional
        Number of hydrogen bonds.
    """

    seq_id: str
    lig_id: str
    pocket_rank: int
    vina_affinity: float
    cnn_pose_score: float
    cnn_affinity: float
    composite_rank_score: float
    overall_rank: int
    n_contacts: Optional[int] = None
    n_hbonds: Optional[int] = None

    def to_dict(self) -> Dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict) -> "RankedPair":
        return cls(**d)
