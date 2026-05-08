# Nature Portfolio Reporting Summary — DAD

> Mr_Prof, 2026-05-08. Drop-in for the official Nature Portfolio Reporting Summary PDF form. Corresponding author transcribes at submission.

## Statistics
- Inferential tests: not used.
- Descriptive statistics: counts, medians, ranges, Pearson r, Spearman rho, AUROC, AUPRC, EF@K.
- Sample sizes: Tier 1 replay n = 6 pairs; RCSB seed n = 16 pairs (4 families).
- Replication: every quantitative claim is re-runnable from 06_Report/Mr_Repro/benchmark/run_replay.py and run_rcsb_redocking.py.
- Randomisation / blinding: not applicable; seed = 0 is fixed for GNINA reproducibility.

## Software
ColabFold >= 1.6.1, GNINA 1.3.2 (master:f23dd2b, built 2025-07-08), P2Rank >= 2.4, DeepTMHMM v1.0.24, Phobius 1.01 (opt-in), RDKit 2024.09+, Open Babel 3.1+, PLIP 2.3+, MMseqs2 release_15, Snakemake >= 8.0, Biopython, ChimeraX >= 1.7, py3Dmol. Pinned CUDA / cuDNN versions in 06_Report/Mr_Repro/tools/gnina_runtime/MANIFEST.md; verified WSL2 + Ubuntu 24.04 runtime in 06_Report/Mr_Pipeline/runtime_freeze.md.

## Materials and reagents
DAD is computational; no physical materials or reagents. Inputs are FASTA + SMILES files and optional precomputed receptor PDB / P2Rank pocket CSVs.

## Antibodies / cell lines / animals / human participants / clinical data / dual-use
Not applicable.

## ChIP-seq / flow cytometry / MRI
Not applicable.

## Data availability
See Publication/submission_metadata/data_availability.md.

## Code availability
See Publication/submission_metadata/code_availability.md.

## Reference-data restrictions
HMDB-derived ligand SMILES (if used) are CC BY-NC 4.0 — non-commercial use only; DAD does not redistribute HMDB data. All other reference databases allow commercial use. PDBbind and CASF-2016 are excluded under the project's free-data policy.

## AI assistance
Software development and manuscript drafting were assisted by large-language-model agents (Anthropic Claude) under human supervision. External quality review by an independent codex agent across six rounds (2026-05-07 to 2026-05-08). All scientific decisions and final manuscript content are the responsibility of the listed human authors. Review history in codex-response/.

## Sampling
- Tier 1 replay: not sampled (all six pairs from supporting paper).
- RCSB seed: curated by family per criteria in Publication/supplementary/SI_text.md SI-M.3.

## Reproducibility footprint
- Tier 1 replay: < 30 seconds, no GPU. 6 / 6 PASS, delta = 0 expected.
- RCSB seed: ~ 4 minutes on RTX 2060. 16 / 16 GNINA execution; 12 / 16 top-pose PASS; 15 / 16 best-of-9 PASS expected.
- All artifact paths and SHA-256 checksums in the corresponding manifest.json files.

## Conditional disclosures (Codex Round 6 unsafe-claim ban list compliance)
- Supporting primary paper status: peer-reviewed and published in *iScience* 2026 (Sung J-Y et al.; PII S2589-0042(26)01130-2; DOI 10.1016/j.isci.2026.<article_id>; URL https://www.cell.com/iscience/fulltext/S2589-0042(26)01130-2). Nature Protocols' supporting-research-paper requirement is satisfied.
- Tier 2 ~ 60-pair generalisation: not completed (planned Phase D / companion-study scope).
- Decoy / enrichment evaluation: not completed (enrichment_metrics.py implemented; dataset run is Phase D).
- Live HTML / .cxc pose visualisation: not produced (Mr_Dock + DAD_protocol.ipynb §11 deliverable).
- Best-of-9 RMSD recovery (15 / 16): search-recovery diagnostic only, not a prospective operational metric. Top-pose (12 / 16) is the operational metric reported.
- Score reproducibility across runtimes: comparable scores within tolerance (Tier 1 acceptance delta_vina <= 0.5, delta_cnn_pose <= 0.05, delta_cnn_aff <= 0.5); exact bitwise identity across GPU architectures is not assumed.
