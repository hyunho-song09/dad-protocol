# DAD — Supplementary Information

> Companion to `Publication/manuscript/DAD_manuscript_v1.md`.
> Mr_Prof, 2026-05-08.
>
> **Drafting rule (same as main)**: every subsection opens with a one-sentence value statement and follows with the quantitative evidence; no descriptive filler. Codex Round 6 unsafe-claim ban list is enforced.

---

## SI-Methods — Detailed methodology — *target ≈ 2000 words*

### SI-M.1 Stage 3 triage rules in full

The Stage 3 triage layer reduces the candidate ORF pool to dock-eligible sequences before any GPU-bound work. Five rules are evaluated in order (`Mr_Bio/rationale.md` §3):

**Rule R1 — Length filter.** Reject ORFs shorter than 50 amino acids. AlphaFold pLDDT becomes meaningful only above ~ 50 residues (Jumper et al. 2021), so shorter sequences cannot produce a reliable structure for downstream pocket detection. Edge case: genuine short binding proteins (lipid-transfer proteins, ribosomal stress peptides) are excluded by default; users can lower the cutoff with `--min-dock-length`.

**Rule R2 — Signal peptide clipping.** If the topology predictor (DeepTMHMM by default; Phobius opt-in) reports a signal peptide, clip the predicted cleavage site and dock the mature form. Signal peptides are removed co-translationally; docking against the unprocessed sequence yields false binding sites on regions that do not exist in the mature protein (Nielsen 2017).

**Rule R3 — Topology-class branching by transmembrane-helix count `nTM`.**

| `nTM` | Class | Action |
|-------|-------|--------|
| 0 | Soluble cytoplasmic / periplasmic | Pass full sequence (path B) |
| 1–2 | Single- or double-pass receptor | Extract extracellular / periplasmic domain (path C) |
| 3–6 | Multi-pass | Extract longest non-TM loop ≥ 60 aa; else FLAG (path D) |
| ≥ 7 | Polytopic integral-membrane | EXCLUDE (path E) |

The MCP case in Tier 1 (NA23_RS01195) exercises path C: the periplasmic sensory domain between TM1 and TM2 is extracted (~ 150 aa, well above the 60 aa minimum) and used as the dock target. Polytopic transporters with `nTM ≥ 7` (e.g., 12-TM major facilitator superfamily) are excluded by default because their extracellular loops are short and disordered, and their ligand-binding faces are lipid-embedded; users can override with `--force-include` for cryptic cytoplasmic sites.

**Rule R4 — Dock-eligible region length.** After R2 / R3 trimming, the remaining sequence must be ≥ 50 aa. Below this, P2Rank pocket detection fails reliably (Krivák & Hoksza 2018).

**Rule R5 — Functional-class boost.** ORFs annotated as known sensors / SBPs / regulators (MCP, CRP/FNR, periplasmic SBP, response regulators) receive a `+1` priority tier. ORFs annotated as DUF / hypothetical are downranked to `verdict = downrank` (still docked but flagged). Annotation source: HMMER vs Pfam by default; Foldseek vs AFDB / PDB optional.

### SI-M.2 External I/O contract for triage

The triage output `triage_report.tsv` carries two columns relevant to downstream stages:

- **`triage_status`** (internal biological label): `PASS / PASS_CLIPPED / FLAG / EXCLUDE`.
- **`verdict`** (external pipeline contract; the only field consumed by Stage 4): `accept / downrank / exclude`.

The mapping is fixed: `PASS` → `accept`; `PASS_CLIPPED` → `accept`; `FLAG` → `downrank`; `EXCLUDE` → `exclude`. This split was introduced in Codex Round 2 to prevent biological labels from accidentally gating downstream stages; the user-facing verdict is now stable across triage tool changes.

### SI-M.3 RCSB co-crystal seed selection

The 16 RCSB co-crystal cases (`build_rcsb_seed_validation.py`) were selected to span the four protein families most relevant to the protocol's user persona (MCP-LBD, CRP/FNR cAMP receptor, periplasmic sugar-binding protein, amino-acid binding protein), with the following filters:

- All four families represented; minimum 2 cases per family where available.
- Resolution ≤ 2.5 Å X-ray (range observed: 0.92 – 2.5 Å; median ~ 1.9 Å).
- Crystallographic ligand chemically valid (atom completeness check; this is what motivated the Round 4 fix replacing 1MDP with 3MBP).
- Receptor PDB CC0; ligand SDF extractable from deposited co-crystal.
- Box centre derived from crystallographic ligand centroid; box size fixed at 22 Å (the auto-box rule's minimum).

The 16 cases are: CRP_SP1_4R8H, CRP_CMP_1I5Z, CRP_CMP_1G6N, CRP_CMP_2CGP (CRP/FNR family); TAR_ASP_2LIG, TAR_ASP_1VLT (MCP-LBD); LAOBP_LYS_1LST, LBP_PHE_1USI, LBP_LEU_1USK (amino-acid binding); GGBP_GLC_2FVY, MBP_GLC_1ANF, MBP_GLC_3MBP, RBSB_RIB_3KSM, GGBP_GLC_2B3B, RBSB_RIB_1DBP, RBSB_RIB_1DRJ (periplasmic sugar-binding).

### SI-M.4 RMSD evaluator — receptor-frame in-place, no Kabsch

The pose-quality metric is receptor-frame in-place heavy-atom RMSD with RDKit symmetry correction. **Kabsch alignment is deliberately not applied.** Kabsch superposition rotates the docked pose onto the crystal pose and minimises RMSD by translation + rotation, removing any receptor-frame placement error. For docking pose evaluation this hides genuine search failures: a pose that is misplaced relative to the receptor by 12 Å (RBSB_RIB_3KSM top pose) can look acceptable post-Kabsch.

The fixed evaluator (R5 fix; `06_Report/Mr_Repro/benchmark/run_rcsb_redocking.py`) computes:

1. Heavy-atom coordinates of the docked pose, in the receptor frame from GNINA's output SDF.
2. Heavy-atom coordinates of the crystallographic ligand, in the same receptor frame from the deposited SDF.
3. Symmetry-corrected RMSD via `rdkit.Chem.AllChem.CalcRMS` with `prbId = -1, refId = -1, maxMatches = 10000` (sufficient for typical drug-like and dipeptide ligands).
4. If RDKit is unavailable, fall back to atom-order in-place RMSD (no symmetry correction).

The output TSV records `top_pose_rmsd_a`, `best_pose_rmsd_a`, `best_pose_index`, `pose_count`, and `rmsd_method` columns. `rmsd_method` is `rdkit_symmetry` for the verified runs.

### SI-M.5 Score parser — first-pose SDF property reading

GNINA's output SDF carries pose properties as `>  <PropertyName>` tags followed by the numeric value on the next line. The pre-R5 parser captured the literal tag string (`<CNNscore>`) when it should have advanced one line and parsed the number. The fixed parser (`_parse_first_pose_props()` in `run_rcsb_redocking.py`) reads `minimizedAffinity`, `CNNscore`, and `CNNaffinity` from the first pose only, taking the value line that follows each tag. All score columns in `redocking_results.tsv` are produced by this parser.

### SI-M.6 Auto-box rule

The Stage 7 box configuration is `box_size = max(22 Å, ligand_max_dim + 10 Å)`, applied isotropically (`size_x = size_y = size_z`). The 22 Å minimum derives from empirical practice in AFPAP / Smiles2Dock and was exercised on all 16 RCSB seed cases without override. Box centre is the P2Rank pocket centroid (live mode) or the crystallographic ligand centroid (RCSB seed mode).

### SI-M.7 Three-runtime equivalence

Three execution paths are documented (`runtime_freeze.md` §1): (A) WSL2 + Ubuntu 24.04 with the project-local CUDA / cuDNN runtime; (B) Google Colab T4 with apt + pip CUDA bootstrap; (C) Docker / Apptainer for HPC. All three use the same canonical command and the same default parameters. **Comparable scores within tolerance are expected; exact bitwise identity across GPU architectures is not assumed**; the Tier 1 acceptance tolerance is Δ_vina ≤ 0.5 kcal mol⁻¹, Δ_cnn_pose ≤ 0.05, Δ_cnn_aff ≤ 0.5.

---

## SI-Notes — Per-case failure analysis — *target ≈ 1500 words*

The four top-pose failures from §10.2 Table 6 are reported here with the per-case detail from `06_Report/Mr_Repro/external_validation/rcsb_seed/failure_analysis.md` (Mr_Repro, 2026-05-08). Verbatim citation of the operative paragraphs follows; this section in the protocol's submitted SI will reproduce them exactly.

### SI-N.1 TAR_ASP_1VLT — within-case ranking failure

**Protein.** TAR (methyl-accepting chemotaxis protein), aspartate-binding domain. **PDB.** 1VLT (resolution 2.2 Å). **Ligand.** L-aspartate (133 Da, sp³ backbone + carboxyl + amino).

**RMSD profile.** Top-pose 3.009 Å (FAIL); best-pose 1.627 Å at index 4 (PASS); 9 generated modes.

**Interpretation.** Search recovered a valid conformation (best RMSD < 2 Å), but CNN scoring placed an alternative orientation first with confidence 0.9549. Aspartate is a small zwitterion in a narrow electrostatic channel; rotational ambiguity of the carboxyl group within the pocket creates multiple chemically equivalent orientations of nearly identical electrostatic energy. The CNN scoring head, trained on CrossDocked2020, has not learned the experimental orientation as the dominant mode for this case.

### SI-N.2 MBP_GLC_1ANF — true search failure

**Protein.** Maltose-binding protein, glucose-binding domain. **PDB.** 1ANF (resolution 1.67 Å). **Ligand.** glucose (180 Da, cyclic polyol).

**RMSD profile.** Top-pose 5.124 Å (FAIL); best-pose 3.464 Å at index 8 (also FAIL — > 2.0 Å threshold); 9 generated modes.

**Interpretation.** No pose in the 9-mode ensemble cleared the 2.0 Å threshold. This is a search failure rather than a ranking issue. MBP undergoes a large Venus-flytrap domain closure on substrate binding; the rigid-receptor docking against the deposited (closed) PDB does not sample the closed-form geometry that the actual binding mode requires. Rigid-receptor GNINA cannot, in principle, rescue this case without explicit receptor-flexibility sampling.

### SI-N.3 MBP_GLC_3MBP — within-case ranking failure (same family as SI-N.2)

**Protein.** Maltose-binding protein. **PDB.** 3MBP (resolution 1.7 Å). **Ligand.** glucose.

**RMSD profile.** Top-pose 5.003 Å (FAIL); best-pose 0.626 Å at index 7 (PASS); 9 generated modes.

**Interpretation.** Search recovered an excellent pose at index 7 (RMSD 0.626 Å), but CNN scoring ranked a 5.003 Å pose first. The contrast with SI-N.2 (same protein family, same ligand, different crystal form) demonstrates that **search and ranking are independent failure modes within a single family**. The protocol's three-score reporting plus the best-of-N column is the only diagnostic that catches both modes; the operational fix (`select_consensus_pose()`) is Phase D scope.

### SI-N.4 RBSB_RIB_3KSM — severe within-case ranking failure

**Protein.** RbsB (ribose-binding protein), substrate-binding domain. **PDB.** 3KSM (resolution 1.9 Å). **Ligand.** β-D-ribose (150 Da, pentose).

**RMSD profile.** Top-pose **12.778 Å** (severe FAIL); best-pose 1.572 Å at index 4 (PASS); 9 generated modes. CNN top-pose confidence 0.8943; Vina score -3.706 kcal mol⁻¹ (the weakest in the 16-case set).

**Interpretation.** The pose ensemble contained a near-native solution (RMSD 1.572 Å), but the top rank is globally misplaced (RMSD > 12 Å, likely a flipped or reversed orientation). Both Vina and CNN top scores are the weakest in the benchmark, indicating low confidence overall. **This is the single case that justifies the §SI-M.4 demand for receptor-frame in-place RMSD**: under Kabsch alignment, the misplaced top pose would have looked acceptable.

### SI-N.5 Family-level pattern and three hypotheses

The failures cluster on periplasmic sugar-binding proteins (4 / 7 family pass-rate). Three non-exclusive hypotheses, drawn from `failure_analysis.md` §"Family-Level Pattern":

1. **Multiple binding modes.** Sugars (glucose, ribose) carry several chemically equivalent hydroxyl groups that can engage symmetric hydrogen-bond networks with the pocket. The energy landscape for these orientations is flatter than for amino-acid or nucleotide ligands; CNN scoring, even when trained on a large set, cannot reliably prefer the experimental orientation when alternates are equally plausible.
2. **Induced-fit / receptor flexibility.** MBP, GGBP, RbsB undergo substantial domain closure and hinge motion on substrate binding. Rigid-receptor docking against a single crystallographic frame cannot model this, systematically excluding the actual bound geometry. The MBP_GLC_1ANF search failure is the dominant manifestation.
3. **CNN training-distribution bias.** GNINA's CNN was trained on CrossDocked2020, which is dominated by rigid-pocket, drug-like-ligand complexes. Periplasmic sugar-binding proteins may be under-represented in the training distribution, making them out-of-distribution for the scoring head and amplifying the within-case ranking failures observed for MBP and RbsB.

### SI-N.6 Operational consequence (cross-references main §10.1, §11.1)

DAD users targeting periplasmic sugar-binding proteins should treat the top-1 pose as advisory and inspect the best-of-N column as a diagnostic — not as a substitute. Concrete mitigations beyond best-of-N reporting (flexible-receptor docking; ensemble docking against MD-derived states; alternative scoring functions; the in-tree `select_consensus_pose()` consensus rule under construction in Phase D) are out of scope for this manuscript.

---

## SI-License — Reference data license audit — *target ≈ 500 words*

DAD's reference-data dependency tree consists entirely of free, redistribution-permissive sources. The only NC-flagged optional source is HMDB; PDBbind and CASF-2016 are deliberately excluded under the project's free-data policy.

### Table SI-L. Reference databases used or planned by DAD.

| Database | License | Commercial use | Redistribution | Citation duty | DAD usage |
|----------|---------|---------------|---------------|---------------|-----------|
| RCSB PDB | CC0 | Allowed | Allowed | Recommended | Tier 1 + RCSB seed structures |
| AlphaFold DB | CC BY 4.0 | Allowed | Allowed | Required | Tier-2 structure augmentation (planned) |
| BindingDB | CC BY 3.0 | Allowed | Allowed (attribution) | Required | Tier-2 affinity (planned) |
| BioLiP2 | Academic free + commercial-permitted (verify) | User-verified | Allowed (snapshot) | Required | Tier-2 binding sites (planned) |
| CrossDocked2020 | CC0 / MIT | Allowed | Allowed | Recommended | Tier-2 pose-quality subset (planned) |
| ChEMBL | CC BY-SA 4.0 | Allowed | Allowed (share-alike) | Required | Tier-2 transporter targets (planned) |
| MetaboLights | CC BY 4.0 | Allowed | Allowed | Required | Optional metabolite SMILES pool (planned) |
| HMDB | CC BY-NC 4.0 | **Not allowed** | Academic only | Required | Optional metabolite SMILES pool (planned, NC-flagged) |
| PoseBusters | MIT | Allowed | Allowed | Recommended | Pose-validity gate (Phase D) |
| **PDBbind** | (excluded) | — | — | — | **Excluded by user policy** |
| **CASF-2016** | (excluded) | — | — | — | **Excluded by user policy** |

### License compliance procedure

`Mr_Repro/data_manifest.tsv` records, per case: source database, license flag, access date, SHA-256 of the downloaded artifact, and (for HMDB) an explicit NC flag column. Tier-1 worked-example data are user-supplied and inherit the user's institutional policy; Tier-2 ingest scripts are scaffolded under `06_Report/Mr_Repro/data/ingest_phase_d/` and are runnable without registration.

The Reporting Summary states that HMDB-derived ligand SMILES, if used, require non-commercial use under CC BY-NC 4.0, and that the protocol does not redistribute HMDB data — the ingest script is a download tool, not a redistribution surface.

### Why PDBbind / CASF were excluded

PDBbind and CASF-2016 are the conventional benchmarks for scoring-function calibration but require institutional-email registration and carry CC BY-NC-SA licenses. Both are direct reproducibility barriers for industry users and for any reader without institutional credentials. The project's free-data policy excludes them; PoseBusters (MIT) supersedes CASF for pose-validity testing in 2024+ literature (Buttenschoen et al. 2024) and is used by DAD in Phase D scope.

---

*End of SI text. Numerical claims in SI-Methods and SI-Notes are reproducible from the corresponding artifacts in `06_Report/`. SI-License Table SI-L is sourced from `06_Report/Mr_Bio/dataset_strategy.md` v2.*
