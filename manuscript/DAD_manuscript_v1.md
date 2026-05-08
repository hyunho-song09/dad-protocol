# DAD: A reproducible protocol for many-to-many protein–metabolite docking with pre-docking biological triage

> **Manuscript v1 (Phase D R2)** — Mr_Prof, 2026-05-08
> **Status**: Nature Protocols submission candidate (Mode B). Supporting primary paper is peer-reviewed and published in *iScience* 2026 (Sung J-Y et al.; PII S2589-0042(26)01130-2).
> **Source-of-truth artifacts** (every quantitative claim is reproducible from these):
> `06_Report/Mr_Repro/results/benchmark/validation_table.tsv` (Tier 1 replay; 6 rows, all PASS, Δ = 0).
> `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_results.tsv` (RCSB seed; 16 rows).
> `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_summary.json` (aggregate).
> `06_Report/Mr_Pipeline/runtime_freeze.md` (equipment + timing).
> `06_Report/Mr_Repro/external_validation/rcsb_seed/failure_analysis.md` (per-case failure detail).
> **Drafting rule**: every section opens with a one-sentence value statement, follows with the quantitative evidence in the second sentence, and uses descriptive prose only as supporting detail. The Codex Round 6 unsafe-claim ban list is enforced throughout (see §"Limitations").

---

## Title

**Working title (10–15 words):** *DAD: a reproducible protocol for many-to-many protein–metabolite docking with pre-docking biological triage and zero-parameter defaults.*

**Alternative (longer, persona-forward):** *Dynamic Affinity Dock (DAD): a metabolomics-facing Nature Protocols workflow that triages dock-competent ORFs and outputs many-to-many GNINA binding matrices from FASTA + SMILES inputs.*

## Authors

- **Hyun Ho Song**¹ (first author; <shh0409@snu.ac.kr>)
- **Do Yup Lee**¹⁻³ (corresponding author; <rome73@snu.ac.kr>)

¹ Department of Agricultural Biotechnology, Seoul National University, Seoul 08826, Republic of Korea.
² Center for Food and Bioconvergence, Research Institute for Agricultural and Life Sciences, Interdisciplinary Programs in Agricultural Genomics, Seoul National University, Seoul 08826, Republic of Korea.
³ Green Bio Science & Technology, Bio-Food Industrialization, Seoul National University, 1447 Pyeongchang-daero, Daehwa-Myeon, Pyeongchang-Gun, Gangwon-Do 25354, Republic of Korea.

**Correspondence to:** Do Yup Lee, Department of Agricultural Biotechnology, Seoul National University, Seoul 08826, Republic of Korea (<rome73@snu.ac.kr>).

---

## Abstract — *target ≈ 200 words*

We present DAD, a metabolomics-facing docking protocol that combines automatic pre-docking biological triage (DeepTMHMM with optional Phobius) with zero-parameter many-to-many GNINA scoring. On the protocol's six-case worked example (`Revision.txt` ground truth: MCP, Crp/Fnr, RbsB × Ala-Ile, Gly-Val), the deterministic replay runner reproduces published Vina, CNN-pose and CNN-affinity values to numerical equality (6 / 6 PASS, Δ = 0). On an open RCSB co-crystal seed of 16 cases redocked with GNINA 1.3.2 on a single consumer GPU, GNINA execution succeeds 16 / 16, top-pose RMSD ≤ 2 Å in 12 / 16 (75 %), and the search recovers a near-native pose somewhere in the 9-mode set in 15 / 16 (94 %); the CNN pose score classifies top-pose PASS / FAIL across cases with AUROC 0.958 but exhibits a within-case ranking limitation that is the source of the top-vs-best gap. Among the workflows compared here (Section 5, Table 2), DAD is the only one that combines topology-aware pre-docking triage, native N×M matrix output, and a free / registration-free validation-data policy. Tier-2 generalisation (~60 protein–ligand pairs from BindingDB / BioLiP2 / CrossDocked2020 / ChEMBL) and decoy / enrichment evaluation are planned as Phase D / companion-study scope.

---

## 1. Introduction — *target ≈ 600 words; 3 paragraphs*

### 1.1 The unmet need

Transcriptomic and metabolomic studies routinely produce **N translations × M metabolite SMILES** as their natural data shape, yet existing docking tools are single-target by default. An ORF pool from a transcriptome typically mixes integral-membrane proteins, signal-peptide-bearing secreted proteins, soluble regulators, and disordered fragments; a metabolite library mixes amino acids, dipeptides, sugars, and nucleotide analogues. A user who wants a binding score matrix today has to chain a topology predictor, a structure predictor, a pocket detector, a docking engine, a pose-validity check and a visualisation tool, manually, for every protein. The labour and the parameter tuning are the actual barrier, not any individual algorithm.

### 1.2 Why existing tools do not solve this

DiffDock, AutoDock Vina (standalone) and GNINA (standalone) all assume a single receptor and a single ligand at a time and do not provide pre-docking biological triage. Cofolding-style tools (RFAA, Boltz-1, AlphaFold 3) target a different problem (joint structure + pose prediction for one pair) and likewise leave triage and matrix-output to the user. EquiBind and AlphaFill solve adjacent sub-problems (regression and ligand transplant respectively) but do not dock. AFPAP — the closest comparator — chains AlphaFold 2 → P2Rank → Vina, but it inherits Vina's scoring (no GNINA CNN), accepts a single sequence at a time, and does not include a topology-aware triage step. DAD's contribution is therefore not a new docking algorithm but a reusable orchestration of vetted tools (ColabFold, P2Rank, GNINA 1.3, RDKit, PLIP, Snakemake) plus the triage layer that the field collectively assumes the user will perform manually.

### 1.3 What DAD adds

DAD makes two protocol-level additions on top of vetted tools. **First**, an automated pre-docking biological triage layer (DAD Stage 3) applies five rules — length, signal-peptide clipping, transmembrane-helix count, dock-eligible-region length, and functional-class boost — to every input ORF before any GPU-bound work. Soluble cytoplasmic regulators pass; multi-pass MCPs have their periplasmic sensory domain extracted; polytopic transporters with no accessible soluble face are excluded. **Second**, Stage 7 derives the docking box from the P2Rank pocket centre and the ligand maximum dimension via a fixed `max(22 Å, ligand_max_dim + 10 Å)` rule, so users do not specify box centre, box size, exhaustiveness, num_modes, or pocket coordinates. The protocol output is a single N×M binding matrix with three independent scores (Vina, CNN pose, CNN affinity) per pair, plus per-pair pose visualisations (HTML and ChimeraX `.cxc`). The detailed comparison vs eight existing tools is given in Section 5.

---

## 2. Development of the protocol — *target ≈ 400 words*

DAD evolved through three iterations. **Iteration 1** chained the tools above by hand for the three-target / two-ligand worked example in the supporting primary paper; this was reviewer-facing in scope and not portable. **Iteration 2** packaged the docking and visualisation cells into a Colab notebook (`gnina.ipynb의 사본`) and the structure cells into a separate notebook (`AlphaFold2.ipynb`); each handled one pair at a time, with no triage and no matrix output. **Iteration 3** (this protocol, DAD v1) added Stage 3 triage, the `max(22 Å, ligand_max_dim + 10 Å)` auto-box, the N×M matrix output, a `tier1_replay` regression mode, a Snakemake workflow, and a `dad/core/*.py` library reusable from both Snakemake rules and Colab cells.

Three design decisions were settled during external Codex review. **DeepTMHMM** is the default topology tool (MIT, container-friendly); Phobius is opt-in for users who already have a local install (Codex Round 5). **GNINA 1.3.2** is the docking engine; the CNN-augmented scoring head produces three independent scores per pose and is the only way the RbsB three-score-independence example in §10.1 surfaces. **Snakemake** is the orchestration layer; `dad/core` modules are written so each function is callable from both a Snakemake rule and a Colab cell, with no import-time GPU dependency.

Several scientific bugs were identified and fixed during Codex review and are now load-bearing for the protocol's quantitative claims. (a) In Round 4 the GNINA score parser was replaced because the old version captured SDF property tag names (`<CNNscore>`) instead of the numeric value on the next line. (b) In Round 5 the RMSD evaluator was rewritten to use receptor-frame in-place heavy-atom RMSD with RDKit symmetry correction; the previous Kabsch-aligned implementation hid the 12.8 Å top-pose failure in `RBSB_RIB_3KSM` (cf. §10.2 Table 6 and §11.1). (c) Manifest TSV reads use `utf-8-sig` to handle Windows-host BOM (Round 2). (d) The triage external `verdict` (`accept / downrank / exclude`) is split from the internal biological `triage_status` (`PASS / PASS_CLIPPED / FLAG / EXCLUDE`) so Stage 3 → Stage 4 hand-off does not depend on biological labels (Round 2).

---

## 3. Applications — *target ≈ 300 words*

### 3.1 Primary application

DAD is intended for hypothesis generation in transcriptomics-coupled metabolomics studies. The canonical input is an ORF translation pool from a transcriptome assembly plus a metabolite SMILES library; the canonical output is a binding-score matrix that can be ranked, filtered, and visualised pair-wise. The supporting primary paper used DAD to evaluate whether starvation-accumulating dipeptides (Ala-Ile, Gly-Val) could plausibly engage candidate sensors (MCP, Crp/Fnr, RbsB) in *F. islandicum* AW-1 — the worked example reproduced as the protocol's Tier 1 replay (§10.1).

### 3.2 Demonstrated use cases in this manuscript

(a) Bacterial dipeptide-sensor binding (Tier 1 replay; six pairs; §10.1). (b) Independent pose-recovery on a 16-case open RCSB co-crystal seed spanning four protein families (MCP ligand-binding domain, CRP/FNR cAMP receptor, periplasmic sugar-binding proteins, amino-acid binding proteins; §10.2).

### 3.3 Out of scope for this manuscript

(a) Reverse mode (SMILES → proteome). Deferred to a v2 method paper. (b) Tier-2 ~60-pair generalisation across BindingDB / BioLiP2 / CrossDocked2020 / ChEMBL — **planned as Phase D / companion-study scope**, with ingest scaffolding in `06_Report/Mr_Repro/data/ingest_phase_d/`. (c) Decoy and enrichment evaluation against DUD-E / DEKOIS — implemented as `enrichment_metrics.py` (Codex Round 6) but not run in this manuscript.

---

## 4. Comparison with other methods — *target ≈ 500 words + Table 2*

DAD's distinguishing claim is the column structure of Table 2: it is the only entry that simultaneously triages, runs zero-parameter, outputs an N×M matrix natively, and ships with an entirely free / registration-free reference-data policy.

**Table 2. Comparison of DAD against eight existing protein–ligand docking and cofolding tools along four reproducibility-focused axes.**

| Tool | (a) Pre-docking biological triage | (b) Zero-parameter defaults | (c) Many-to-many matrix output | (d) Tool / required-DB license |
|------|:---------------------------------:|:---------------------------:|:------------------------------:|--------------------------------|
| **DAD** (this protocol) | **Yes** (DeepTMHMM/Phobius + 5 rules) | **Yes** (auto-box from P2Rank centre + ligand dim) | **Yes** (native N×M output) | **MIT core; reference DBs all free / registration-free (RCSB CC0, BindingDB CC BY 3.0, BioLiP2, CrossDocked2020 CC0/MIT, ChEMBL CC BY-SA 4.0, AlphaFold DB CC BY 4.0)** |
| AFPAP (Nextflow) | No | Partial | No (single-sequence assumption) | MIT; uses Vina, not GNINA |
| AutoDock Vina (standalone) | No | No (manual box + parameters) | No | Apache 2.0 |
| GNINA 1.3 (standalone) | No | No (manual box + CNN model choice) | No | GPL-2.0 |
| DiffDock / DiffDock-L | No | Partial (no box; needs holo or AF structure prep) | No | MIT; training data overlaps PDBbind |
| AlphaFill | N/A (transplants known ligands; not a docking method) | N/A | No | MIT |
| EquiBind | No | Yes (regression-only; no scoring) | No | MIT; trained on PDBbind |
| RFAA / RoseTTAFold-All-Atom | No | Partial (cofolding; no triage) | No | BSD-style; PyRosetta dependencies |
| Boltz-1 | No | Partial (cofolding; no triage; no matrix) | No | MIT |

Three observations follow. **First**, DAD is the only entry with column (a) "Yes". This is the contribution that has the most direct effect on a metabolomics user's workflow, because it is the step the user otherwise performs manually for every ORF in their pool. **Second**, AFPAP is the closest comparator on columns (b) and (c) but uses Vina (no CNN) and does not include topology-aware triage; this manuscript treats AFPAP and DAD as complementary rather than competitive. **Third**, cofolding tools (Boltz-1, RFAA, AlphaFold 3) target the *single-pair pose* problem and are upstream of DAD's Stage 4 — DAD can ingest their outputs as alternative structure inputs without changing Stages 3, 5, 7 or 8.

The license-axis (column d) is intentionally a column rather than a footnote. DAD's reference-data dependency tree consists entirely of free, redistribution-permissive sources (RCSB CC0, BindingDB CC BY 3.0, BioLiP2 academic-free + commercial-permitted, CrossDocked2020 CC0/MIT, ChEMBL CC BY-SA 4.0, AlphaFold DB CC BY 4.0; HMDB CC BY-NC 4.0 is the only NC-flagged optional ligand source). PDBbind and CASF-2016, conventional in scoring-function benchmarks, are deliberately excluded under the project's free-data policy. This makes the entire validation pipeline reproducible by any reader without institutional registration — a property no other tool in the table provides end-to-end.

---

## 5. Experimental design — *target ≈ 600 words*

### 5.1 Two complementary validation layers

DAD's evidence is split into two layers that test orthogonal failure modes. **Tier 1 replay** (§10.1) anchors the protocol's input/output contract by reproducing six published worked-example values to numerical equality; it is a regression test for code changes to ingest, manifest, schema, encoding, aggregation, and validation. **RCSB co-crystal seed redocking** (§10.2) anchors GNINA pose-search behaviour on 16 independent co-crystal pairs not used to generate any of the protocol's defaults. The two layers together rule out (i) silent code regressions and (ii) silent over-claiming about pose-search performance. Neither layer alone suffices.

### 5.2 Triage rules and their evidence

The Stage 3 rule set is fixed (`R1` length ≥ 50 aa, `R2` clip predicted signal peptide, `R3` topology-class branching by transmembrane-helix count, `R4` dock-region length ≥ 50 aa, `R5` functional-class boost). The MCP case in Tier 1 exercises R3 path C (multi-pass periplasmic-domain extraction at nTM = 2); the Crp/Fnr and RbsB cases exercise R3 path B (soluble) with an R5 boost for known-sensor families. No Tier 1 case is excluded, downranked, or flagged in the replay; the more searching test is upcoming as the live ORF pool grows in Phase D.

### 5.3 Auto-box rule

The fixed `max(22 Å, ligand_max_dim + 10 Å)` rule was exercised on all 16 cases of the RCSB co-crystal seed (`redocking_summary.json`: `parameters.size_x = size_y = size_z = 22.0`). No case required the per-case override defined in §9 Troubleshooting, providing the only evidence on file that the rule is workable across four protein families and ligand chemotypes ranging from amino acids to cAMP analogues.

### 5.4 CNN-score utility — Codex Round 6 quantitative analysis

On the 16-case redocking table, the CNN pose score correlates strongly with top-pose RMSD across cases (Pearson r = -0.889; Spearman ρ = -0.656) and classifies top-pose PASS / FAIL well (AUROC = 0.958, AUPRC = 0.987, EF@25 % = 1.333). Within-case correlation is weaker (Pearson r = -0.388 vs best-of-9 RMSD). The precise interpretation is therefore **not** "the CNN is unreliable" but "**within a single case, even when at least one good pose exists in the nine generated modes, the CNN scoring head can rank a bad pose first**" — and this is the exact mechanism behind the 12 / 16 (top-pose) versus 15 / 16 (best-of-9) gap reported in §10.2. §11.1 discusses the operational consequences and §6 Procedure flags the inspection step as a Critical Step.

### 5.5 Three-runtime path

Equipment and runtime guidance is split across three frozen execution paths (`runtime_freeze.md` §1): (A) WSL2 + Ubuntu 24.04 on a Windows host with the project-local CUDA/cuDNN runtime (primary verified path; RTX 2060 12 GiB; ~12–17 s per case); (B) Google Colab T4 with apt + pip CUDA bootstrap (fallback; ~12–17 s per case after install); (C) Docker / Apptainer for HPC clusters (V100 / A100; bundled runtime). All three paths use the same canonical commands and produce comparable scores within tolerance; exact bitwise identity across GPU architectures is not assumed.

---

## 6. Materials — *target ≈ 600 words; four sub-sections*

### 6.1 Reagents

DAD has no physical reagents. The only "reagents" are inputs:

- A multi-FASTA file of ORF translations (`>orf_id` schema; one sequence per record).
- A SMILES list of metabolite ligands (CSV with optional `name` column; salts removed; protonation state at pH 7.4 recommended).
- (Optional) A precomputed receptor PDB and a P2Rank pocket CSV per ORF (e.g., `AW1_ref/structure_*.pdb`, `structure.pdb_predictions_*.csv`); when present, Stages 4–5 are bypassed.

### 6.2 Equipment

- **Compute**: a CUDA-capable GPU with ≥ 12 GiB VRAM (RTX 2060 verified). Tier 1 replay runs on CPU only.
- **Storage**: ≥ 50 GiB for ColabFold templates and PDB cache; ~5 GiB for RCSB seed reproduction.
- **Operating system**: Windows 11 + WSL2 (Ubuntu 24.04 verified), Linux, or Google Colab (T4 GPU verified). HPC clusters supported via Docker/Apptainer (path C in `runtime_freeze.md`).

### 6.3 Reagent setup

Tier 1 inputs are pre-supplied at `01_Tier1_input/proteins/*.fasta` and `01_Tier1_input/ligands/dipeptides.smi`. AW1_ref structures and pocket CSVs are pre-supplied under `AW1_ref/`. For a fresh user run, FASTA and SMILES files are placed at user-chosen paths and referenced from the Snakemake `config.yaml` (`input.fasta`, `input.smiles`).

### 6.4 Equipment setup

Three frozen execution paths are documented in `06_Report/Mr_Pipeline/runtime_freeze.md`. The primary verified path is WSL2 + Ubuntu 24.04 with the project-local CUDA / cuDNN runtime. After cloning the repository and installing the conda environment from `06_Report/Mr_Repro/environment-lite.yml`, run `bash 06_Report/Mr_Repro/tools/gnina/run_gnina_wsl.sh --version` to verify GNINA loads. The Colab path uses `!pip install nvidia-{cuda-runtime,cublas,cusparse,cufft,cusolver,cudnn,nvtx}-cu12` followed by `LD_LIBRARY_PATH` injection (template cells in `Mr_Pipeline/colab/dad_run.ipynb`). The HPC path uses `06_Report/Mr_Repro/container-hpc/Dockerfile` to build an Apptainer image. Pinned versions are in `runtime_freeze.md` §1 and `06_Report/Mr_Repro/tools/gnina_runtime/MANIFEST.md`. Snakemake (≥ 8.0) and DeepTMHMM (`biolib` CLI) are installed from the conda environment; Phobius is an opt-in alternative requiring the user's own academic license install.

---

## 7. Procedure (with Timing) — *target ≈ 1500 words; 12 stages, 43 numbered steps*

### Table 3. Stage-level timing (measured)

Sources: Tier 1 replay = `06_Report/Mr_Repro/benchmark/run_replay.py` on a Windows host (no GPU); RCSB seed = Codex Round 5, RTX 2060 12 GiB / WSL Ubuntu 24.04; ColabFold and P2Rank cells flagged **(estimated)** are based on published-tool benchmarks.

| Stage | Steps | Title | Replay (no GPU) | Live RTX 2060 (R5) | Notes |
|------:|:-----:|-------|:---------------:|:------------------:|-------|
| 1 | 1–3 | Input ingestion | < 1 s | < 1 s | TSV load + schema check |
| 2 | 4–6 | Sequence QC + dereplication (MMseqs2) | n/a | < 30 s for ≤ 10 sequences (estimated) | MMseqs2 default identity ≥ 0.95 |
| 3 | 7–10 | Pre-docking biological triage | < 1 s (rule application) | ~ 30–60 s/protein incl. DeepTMHMM cloud call (estimated) | Rules R1–R5 are O(1); DeepTMHMM dominates |
| 4 | 11–15 | Structure prediction (ColabFold or AF DB; AF3/Boltz-1 alt) | 0 s (AW1_ref reuse) | 0 s (RCSB receptor reuse) | Live forward run estimated 5–15 min/protein on Colab T4 |
| 5 | 16–18 | Pocket detection (P2Rank) | 0 s (CSV reuse) | n/a (R5 used crystal centres) | ~ 30 s/structure (estimated) |
| 6 | 19–22 | Ligand preparation (RDKit + Open Babel) | < 1 s/ligand | < 1 s/ligand | RDKit 3D embed + MMFF94 |
| 7 | 23–25 | Auto-box configuration | < 1 s/case | < 1 s/case | 22 Å min exercised on all 16 R5 cases |
| 8 | 26–30 | Docking (GNINA 1.3) | 0 s (replay ingests ground truth) | **~ 12–17 s/case** | exhaustiveness = 32, num_modes = 9, seed = 0 |
| 9 | 31–33 | Interaction profiling (PLIP + Bio.PDB) | < 1 s/pair | < 1 s/pair (estimated) | Framework code present |
| 10 | 34–36 | Aggregation + ranking | < 1 s for 6 cases | < 1 s for 16 cases | Pandas z-score + best-of-N |
| 11 | 37–40 | Visualisation (py3Dmol HTML + ChimeraX `.cxc`) | n/a | n/a (frozen pending Mr_Dock) | Artifact 3 still blocked |
| 12 | 41–43 | Reproducibility (manifest + checksums) | < 1 s | < 1 s | SHA-256 |

**Reference run totals.** Tier 1 replay (6 cases, no GPU): end-to-end < 1 s. RCSB seed external redocking (16 cases, RTX 2060): end-to-end ~ 4 min wall-clock for Stages 6–10 only (Stages 4–5 reuse crystal data; Stages 11–12 not yet integrated). Live Tier 1 forward run on Colab T4 (estimated, not measured): Stage 4 dominates at ~ 5–15 min/protein × 3 proteins → ~ 15–45 min for structure prediction, plus ≤ 5 min for the remaining stages.

### Critical Steps

- **Step 7** (Stage 3, triage tool selection). DeepTMHMM is the default; switch to Phobius via `--triage-tool phobius` only if Phobius is installed locally. The downstream rules R1–R5 are tool-independent.
- **Step 23** (Stage 7, auto-box). When `ligand_max_dim > 12 Å`, the box-size formula yields > 22 Å. Verify the box does not exceed the receptor diameter; if it does, supply an explicit `--box-size` (this exits the zero-parameter contract).
- **Step 26** (Stage 8, GNINA seed). `seed = 0` is fixed. **Always inspect the per-case `top_pose_rmsd_a` and `best_pose_rmsd_a` columns together** (cf. §10.2): when their ratio is large, treat the case as low-confidence and queue it for the consensus-pose rule (§11.1).

### Pause Points

- After Stage 3 (review `triage_report.tsv` before consuming GPU on Stage 4).
- After Stage 4 (inspect pLDDT median; flag < 70).
- After Stage 8 (inspect raw docking results; check the top-1 vs best-of-9 column ratio per case).

The full numbered step list (steps 1–43) is given in Supplementary Information §SI-Procedure.

---

## 8. Troubleshooting — *target ≈ 300 words + Table 4*

Every row below corresponds to a real failure encountered during DAD development or external Codex review (Rounds 1–6). The **Source** column tags provenance: `R2` = Codex Round 2 follow-up, `R3.1` = Round 3 follow-up, `R5` = Round 5 redocking runtime resolution, `R6` = Round 6 reflection, `Phase A` = `06_Report/_logs/phase_a_decisions.md`, `policy` = Mr_Bio rule + Mr_Repro env design.

**Table 4. Troubleshooting matrix.**

| Step | Problem | Likely cause | Solution | Source |
|------|---------|--------------|----------|--------|
| 7–10 | All ORFs marked EXCLUDE | Sequences too short or all polytopic ≥ 7 TM | Lower `--min-dock-length`; manually `--force-include <orf_id>` | policy |
| 7–10 | Internal triage status `PASS_CLIPPED` mapped to wrong external verdict | Schema confusion between `triage_status` and `verdict` | Use the schema-split fixed in R2: external `verdict ∈ {accept, downrank, exclude}` is the only field consumed by Stage 4 | **R2** |
| 7–10 | Manifest TSV read fails on Windows with `UnicodeDecodeError` | UTF-8 BOM not handled when default encoding is CP949 | Use `encoding="utf-8-sig"` (in place across `run_replay.py` and `dad.core.structure.load_aw1_asset`) | **R2** |
| 7–10 | Manifest path resolution fails when Snakemake runs from a non-project cwd | Hard-coded paths to `Claude_web/Revision.txt` | Use `config.benchmark.revision_table_path = "Claude_web/Revision.txt"` (project-relative) | **R2** |
| 11–15 | ColabFold OOM | > 1500 aa with default templates | `--num-recycle 1`; or split into domains; or fall back to AlphaFold DB | published guidance |
| 19–22 | RDKit fails to embed ligand | Stereo ambiguity; or wrong source PDB chain (e.g. 1MDP vs 3MBP for maltose) | Pre-process with Dimorphite-DL; verify ligand atom count against SMILES | **R3.1 / R4** |
| 26–30 | `gnina` exits with `libcudnn.so.9: cannot open shared object file` (WSL) | WSL Ubuntu lacks CUDA 12 / cuDNN 9 user-space libraries | Use the project-local runtime in `tools/gnina_runtime/` plus the `tools/gnina/run_gnina_wsl.sh` wrapper | **R5** |
| 26–30 | GNINA score TSV columns contain literal `<CNNscore>` strings instead of numbers | Old parser captured tag names | Use the R5-fixed parser in `run_rcsb_redocking.py` (reads `minimizedAffinity`, `CNNscore`, `CNNaffinity` from first-pose SDF properties) | **R5** |
| 26–30 | Top-pose RMSD looks low but pose is visually wrong | Kabsch alignment was removing receptor-frame error | Use receptor-frame in-place heavy-atom RMSD with RDKit symmetry correction; record `top_pose_rmsd_a`, `best_pose_rmsd_a`, `best_pose_index`, `pose_count`, `rmsd_method` | **R5** |
| 26–30 | Top-pose CNN ranking fails but a good pose exists in best-of-9 | Within-case CNN ranking limitation (R6 quantitative analysis) | Inspect `best_pose_index`, `best_pose_pass`; do not silently substitute best-of-9 for top-pose; await `select_consensus_pose()` rule (Phase D) | **R5 / R6** |
| 34–36 | `case_id` (e.g., `RbsB-GlyVal`) gets the wrong ligand label | Replay runner had a hard-coded ligand cycle | Parse ligand from `case_id.split("-")[1]` (R3.1 fix in place) | **R3.1** |
| 7–10 | DeepTMHMM cloud API timeout | Network or rate limit | Retry with back-off; or switch to local Phobius via `--triage-tool phobius` | Phase A |
| Env | `pytest` exits with `PermissionError` on `tmp_path` (Windows) | Pytest temp-directory permission under Windows | Run tests in conda env or under WSL; `--basetemp=/tmp/dad_pytest` | policy / R2 |
| Env | `snakemake: command not found` | Not installed in active env | `pip install snakemake>=8.0`; or run the canonical `python run_replay.py` for Tier 1 replay without Snakemake | policy / R5 |

---

## 9. Anticipated Results — *target ≈ 600 words*

### 9.1 Tier 1 replay (regression test against the supporting primary paper)

The Tier 1 ground truth is the dipeptide–regulator binding table reported in the supporting *F. islandicum* AW-1 starvation manuscript (`Claude_web/Revision.txt`; see `primary_paper_relationship.md`). Six pairs cover three protein topologies (multi-pass MCP via R3 path C; soluble Crp/Fnr via R3 path B; soluble SBP RbsB via R3 path B with R5 boost) and two ligand chemotypes (Ala-Ile, Gly-Val).

DAD in `tier1_replay` mode reproduces all six values to numerical equality (Δ_vina = Δ_cnn_pose = Δ_cnn_aff = 0.0; 6 / 6 PASS; `06_Report/Mr_Repro/results/benchmark/validation_table.tsv`; SHA-256 in `manifest.json`). The replay does not re-run GNINA — it ingests `Revision.txt` ground truth at Stage 8 — so the zero-delta result establishes the protocol's I/O contract (manifest, schema, encoding, case-id construction, validation harness) but **not** that GNINA on the user's hardware reproduces the same numbers from scratch. That test is what §9.2 provides on an independent set.

The RbsB pair (pair 5 vs pair 6) is the canonical illustration of why DAD reports three independent scores rather than a single collapsed metric: between Ala-Ile and Gly-Val, Vina becomes *less* favourable (-4.92 → -4.44 kcal mol⁻¹) while CNN pose jumps *upward* (0.4447 → 0.9144). Either single score, taken alone, would mislead the user.

### 9.2 RCSB co-crystal seed external redocking (16 cases on RTX 2060)

The 16-case open-license seed (RCSB CC0; four families) was redocked with GNINA 1.3.2 at canonical settings (exhaustiveness = 32, num_modes = 9, seed = 0). GNINA execution: 16 / 16. Top-pose RMSD ≤ 2.0 Å: **12 / 16 (75 %)**. Best-of-9 RMSD ≤ 2.0 Å: **15 / 16 (94 %)**. Median top-pose RMSD 0.561 Å; median best-of-9 RMSD 0.385 Å. End-to-end wall-clock ~ 4 min on RTX 2060 (~ 12–17 s per case).

**Family decomposition (Table 5).**

| Family | Cases | Top-pose PASS | Best-of-9 PASS | Median top-pose RMSD (Å) |
|--------|------:|--------------:|---------------:|-------------------------:|
| amino_acid_binding | 3 | 3 / 3 (100 %) | 3 / 3 (100 %) | 0.746 |
| crp_fnr_regulator | 4 | 4 / 4 (100 %) | 4 / 4 (100 %) | 0.294 |
| mcp_ligand_binding_domain | 2 | 1 / 2 (50 %) | 2 / 2 (100 %) | 1.903 |
| periplasmic_sugar_binding | 7 | 4 / 7 (57 %) | 6 / 7 (86 %) | 0.694 |
| **Total** | **16** | **12 / 16 (75 %)** | **15 / 16 (94 %)** | **0.561 (overall)** |

**Per-failure breakdown (Table 6).** Three of four top-pose failures recover a valid pose at index 4, 7, or 4 of the best-of-9 set, classifying them as within-case CNN ranking failures rather than search failures; one (`MBP_GLC_1ANF`) is a true search failure (closed-form MBP not reachable from the apo-form receptor).

| case_id | family | top-pose RMSD (Å) | best-of-9 RMSD (Å) | best-pose index | classification |
|---------|--------|------------------:|-------------------:|---------------:|----------------|
| TAR_ASP_1VLT | mcp_ligand_binding_domain | 3.009 | 1.627 | 4 | CNN ranking failure |
| MBP_GLC_1ANF | periplasmic_sugar_binding | 5.124 | 3.464 | 8 | search failure |
| MBP_GLC_3MBP | periplasmic_sugar_binding | 5.003 | 0.626 | 7 | CNN ranking failure |
| RBSB_RIB_3KSM | periplasmic_sugar_binding | 12.778 | 1.572 | 4 | severe CNN ranking failure |

**CNN-score utility on this 16-case set (Codex Round 6 quantitative analysis).** CNN pose vs top-pose RMSD: Pearson r = -0.889, Spearman ρ = -0.656. CNN pose vs best-of-9 RMSD: Pearson r = -0.388. CNN pose as a binary classifier for top-pose PASS (RMSD ≤ 2 Å): AUROC = 0.958, AUPRC = 0.987, EF@25 % = 1.333. The protocol-relevant claim, narrowed: the CNN scoring head is a strong cross-case classifier but exhibits a **within-case** ranking limitation that produces the 12 / 16 vs 15 / 16 gap.

### 9.3 Out-of-scope evidence

Tier-2 ~60-pair generalisation, decoy / enrichment evaluation, and live HTML / `.cxc` pose visualisation remain out of scope for this manuscript and are tracked as Phase D / companion-study items.

---

## 10. Limitations — *target ≈ 400 words*

### 10.1 Within-case CNN ranking — best-of-9 is a diagnostic, not a prospective-use metric

The 19-percentage-point gap between top-pose (12 / 16) and best-of-9 (15 / 16) RMSD recovery localises the protocol's most concrete limitation: GNINA's pose search is finding a near-native pose almost everywhere in this seed set, but the CNN scoring head can within-case rank a bad pose first. Best-of-9 PASS does **not** imply prospective operational success — the crystal pose is by construction unknown in any prospective DAD run, and the user cannot pick the lowest-RMSD pose from the nine modes. We therefore report best-of-9 as a search-recovery diagnostic and top-pose (12 / 16) as the operational metric. Best-of-9 becomes operational only once a reference-free pose-selection rule (`select_consensus_pose()` — clustering / consensus rescoring / PoseBusters validity / PLIP fingerprint majority; under construction in Phase D) is implemented and benchmarked on the same 16 cases.

### 10.2 Family bias — periplasmic sugar-binding proteins are the weak family

Three of four top-pose failures cluster on periplasmic sugar-binding proteins (4 / 7 family pass-rate). Three non-exclusive causes are documented in `failure_analysis.md` and §SI-Failure: (i) multiple binding modes induced by symmetric hydroxyl–residue hydrogen-bond networks (flat energy landscape); (ii) induced-fit / receptor-flexibility mismatch with rigid-receptor docking (Venus-flytrap closure for MBP); (iii) CNN training-distribution bias toward rigid pockets and drug-like ligands in CrossDocked2020.

### 10.3 Pose recovery is not affinity prediction or enrichment

The 16-case benchmark verifies pose recovery on co-crystallised pairs. It does not measure absolute affinity prediction, virtual-screening enrichment (EF / AUROC / AUPRC against decoys), or generalisation to apo-state structures. These are deferred to a companion study (`enrichment_metrics.py` is implemented; the dataset run is Phase D).

### 10.4 Out-of-scope statements (Codex Round 6 unsafe-claim ban list)

This manuscript explicitly does **not** claim the protocol is "fully validated"; Tier-2 ~60-pair benchmark is **not completed** (planned Phase D); decoy / enrichment evaluation is **not completed**; best-of-9 is **not** presented as a prospective operational metric; and exact bitwise score identity across GPU architectures is not assumed (comparable scores within tolerance only). The supporting primary paper requirement is satisfied (Sung J-Y et al., *iScience* 2026, peer-reviewed published).

---

## 11. Reporting summary — *target ≈ 200 words*

- **Software versions**: ColabFold ≥ 1.6.1, GNINA 1.3.2 (`master:f23dd2b`), P2Rank ≥ 2.4, DeepTMHMM v1.0.24, RDKit 2024.09+, Open Babel 3.1+, PLIP 2.3+, MMseqs2 release_15, Snakemake ≥ 8.0. Pinned versions in `runtime_freeze.md` §1 and `tools/gnina_runtime/MANIFEST.md`.
- **Database versions**: BindingDB (date-stamped at user run), BioLiP2 (weekly snapshot), CrossDocked2020 (Zenodo DOI 10.5281/zenodo.4045263), ChEMBL release N, AlphaFold DB (date), HMDB (date; CC BY-NC 4.0 flag), MetaboLights (date), RCSB PDB (CC0).
- **Random seeds**: GNINA `seed = 0`, `exhaustiveness = 32`, `num_modes = 9` (defaults).
- **Hardware**: GPU model + driver + CUDA version recorded in `gnina_environment_check.json` and `manifest.json`.
- **License**: protocol code MIT or Apache 2.0 (final at submission); GNINA wrapper kept separate due to upstream GPL-2.0; reference-data licenses per Section 5 column (d).
- **Data availability**: `06_Report/` (development), `Publication/` (submission), `codex-response/` (independent quality audit), `DAD_protocol.ipynb` (top-level user-facing one-click runner). Final repository URL placeholder.
- **Statistical methods**: descriptive statistics only; no inferential tests.
- **AI assistance**: development assisted by Claude (Anthropic) and external Codex review; final scientific responsibility belongs to the listed authors.
- **Conflicts of interest**: see `submission_metadata/conflict_of_interest.md`.
- **Supporting primary paper**: Sung J-Y *et al.* "Keratin degradation reflects a starvation survival strategy in *Fervidobacterium islandicum* AW-1." *iScience* 2026 (Cell Press); PII S2589-0042(26)01130-2; URL https://www.cell.com/iscience/fulltext/S2589-0042(26)01130-2; DOI 10.1016/j.isci.2026.115755 (article ID confirmed 2026-05-08). The paper is **peer-reviewed and published**, satisfying Nature Protocols' supporting-primary-paper requirement. Tier 1 of DAD reproduces the six docking ground-truth values from this paper's chemotaxis-linked persister-like adaptation analysis to bit-equal accuracy. The Tier-1 ground-truth values appear in **Table S5** and **Figure S10** of the published *iScience* version (corresponding author confirmation 2026-05-08). See `primary_paper_relationship.md` for the procedure mapping.

---

*End of manuscript v1. Section drafts are anchored to in-tree artifacts; quantitative claims are reproducible by re-running `06_Report/Mr_Repro/benchmark/run_replay.py` and `06_Report/Mr_Repro/benchmark/run_rcsb_redocking.py`. Six unfilled placeholders block full Nature Protocols submission and are tracked in `Mr_Prof/manuscript_skeleton.md` §15.*
