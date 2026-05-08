# Data Availability Statement

> Mr_Prof, 2026-05-08. Drop-in for Nature Protocols Reporting Summary §"Data Availability".

---

All data underlying the quantitative claims of this manuscript are publicly available through one of three channels: (i) the DAD project repository, (ii) the supporting primary research paper, or (iii) the public reference databases listed below.

## Project repository (primary location)

The DAD project repository contains all in-tree validation artifacts, runner scripts, and per-case outputs.

- **Repository URL**: https://github.com/hyunho-song09/dad-protocol (corresponding author has completed the GitHub OAuth step; team-lead has proposed `dad-protocol` as a public repository — final URL pending). Proposed git tag at submission: `v0.1.0`.
- **Zenodo archive**: on hold per corresponding author (2026-05-08); to be minted at primary-paper acceptance.
- **Working tree root** (during development): `d:\project\experiment\DAD\` (mirrored to repository at submission).

### Tier-1 replay artifacts

| File | Contents |
|------|----------|
| `06_Report/Mr_Repro/results/benchmark/validation_table.tsv` | 6 rows; expected vs observed Vina, CNN pose, CNN affinity; Δ = 0; PASS / FAIL column |
| `06_Report/Mr_Repro/results/report/docking_master.csv` | aggregated 6-pair scores |
| `06_Report/Mr_Repro/results/manifest.json` | tool versions + SHA-256 checksums |
| `06_Report/Mr_Repro/aw1_ref_manifest.tsv` | provenance for AW1_ref structures and pocket CSVs |
| `Claude_web/Revision.txt` | supporting-paper ground-truth source (verbatim Reviewer #2 response excerpt) |

### RCSB co-crystal seed redocking artifacts

| File | Contents |
|------|----------|
| `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_results.tsv` | 16 rows; per-case Vina, CNN pose, CNN affinity, top-pose RMSD, best-of-9 RMSD, best-pose index |
| `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_summary.json` | aggregate statistics; parameters; rmsd_method `rdkit_symmetry` |
| `06_Report/Mr_Repro/external_validation/rcsb_seed/gnina_environment_check.json` | GNINA binary + CUDA / cuDNN library versions |
| `06_Report/Mr_Repro/external_validation/rcsb_seed/prepared/<case_id>/` | per-case prepared receptor PDB and crystal-ligand SDF |
| `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking/<case_id>/` | per-case docked SDFs and GNINA logs |
| `06_Report/Mr_Repro/external_validation/rcsb_seed/failure_analysis.md` | per-case interpretation (4 failures + 3 family hypotheses) |

### Equipment + runtime

| File | Contents |
|------|----------|
| `06_Report/Mr_Pipeline/runtime_freeze.md` | three frozen runtimes (WSL2 Ubuntu 24.04, Colab T4, Docker / Apptainer) with pinned versions |
| `06_Report/Mr_Repro/tools/gnina_runtime/MANIFEST.md` | CUDA 12 / cuDNN 9 user-space pip-wheel pin list |
| `06_Report/Mr_Repro/tools/gnina/run_gnina_wsl.sh` | WSL wrapper that sets `LD_LIBRARY_PATH` |

### User-facing entry point

`DAD_protocol.ipynb` (top-level) — Colab one-click runner. Re-runs the Tier 1 replay and (optionally) the RCSB seed redocking; emits the same artifacts listed above.

## Independent quality audit

`codex-response/` — six-round Codex external review history (Rounds 1 – 6). Read-only; included for reviewer-facing transparency on the project's bug-fix and wording-correction trail. Notable items: Round 4 score-parser fix, Round 5 RMSD evaluator fix and WSL CUDA / cuDNN runtime resolution, Round 6 within-case ranking quantitative analysis (Pearson r, AUROC, EF@25 %).

## Reference databases (public)

DAD's reference data dependency tree consists entirely of free, redistribution-permissive sources. Per `Publication/supplementary/SI_text.md` Table SI-L:

| Database | URL | License |
|----------|-----|---------|
| RCSB PDB | `https://www.rcsb.org/` | CC0 |
| AlphaFold DB | `https://alphafold.ebi.ac.uk/` | CC BY 4.0 |
| BindingDB | `https://www.bindingdb.org/` | CC BY 3.0 |
| BioLiP2 | `https://zhanggroup.org/BioLiP/` | Academic free + commercial-permitted (verify per use) |
| CrossDocked2020 | Zenodo DOI 10.5281/zenodo.4045263 | CC0 / MIT |
| ChEMBL | `https://www.ebi.ac.uk/chembl/` | CC BY-SA 4.0 |
| MetaboLights | `https://www.ebi.ac.uk/metabolights/` | CC BY 4.0 |
| HMDB | `https://hmdb.ca/` | CC BY-NC 4.0 (NC-flagged in `data_manifest.tsv`) |
| PoseBusters | `https://github.com/maabuu/posebusters` | MIT |

PDBbind and CASF-2016 are **excluded** under the project's free-data policy.

## Supporting primary research paper

The Tier-1 ground-truth values (six dipeptide × regulator pairs) are reproduced from Sung J-Y *et al.*, "Keratin degradation reflects a starvation survival strategy in *Fervidobacterium islandicum* AW-1," *iScience* 2026 (Cell Press); PII S2589-0042(26)01130-2; URL https://www.cell.com/iscience/fulltext/S2589-0042(26)01130-2 (DOI 10.1016/j.isci.2026.<article_id>; exact article id to be confirmed by the corresponding author from the iScience landing-page metadata). The six values appear in the chemotaxis-linked persister-like adaptation analysis of the published version and are reproduced verbatim in `Claude_web/Revision.txt`. The exact figure / supplementary identifier inside the published *iScience* version is being confirmed by the corresponding author and will be substituted in the camera-ready submission. See `06_Report/Mr_Prof/primary_paper_relationship.md` for the full procedure-mapping and citation plan. The supporting primary paper is **peer-reviewed and published**, satisfying Nature Protocols' supporting-research-paper requirement.

## Restrictions

- HMDB-derived ligand SMILES, if used, are CC BY-NC 4.0 — non-commercial use only. The DAD project does not redistribute HMDB data; the ingest script is a download tool only.
- All other reference-data sources allow commercial use.
