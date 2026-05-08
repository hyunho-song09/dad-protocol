# Cover Letter — Nature Protocols Submission

> Mr_Prof, 2026-05-08. Editor-facing draft. ~ 400 words.

---

**To:** The Editors, Nature Protocols
**From:** Do Yup Lee (corresponding author; <rome73@snu.ac.kr>), Department of Agricultural Biotechnology, Seoul National University, Seoul 08826, Republic of Korea, on behalf of: Hyun Ho Song (first author) and Do Yup Lee.
**Date:** [TBD — corresponding author to set at submission]
**Re:** Submission of *DAD: a reproducible protocol for many-to-many protein–metabolite docking with pre-docking biological triage and zero-parameter defaults*

---

Dear Editors,

We are writing to submit DAD (Dynamic Affinity Dock) for consideration at Nature Protocols. DAD is a reproducible end-to-end protocol that allows transcriptomic and metabolomic users to obtain a many-to-many protein–metabolite binding-score matrix from FASTA + SMILES inputs without parameter tuning. It packages the methodology used in our supporting peer-reviewed primary research paper on starvation-induced peptide signalling in *Fervidobacterium islandicum* AW-1 (Sung J-Y *et al.*, *iScience* 2026, Cell Press; PII S2589-0042(26)01130-2; URL https://www.cell.com/iscience/fulltext/S2589-0042(26)01130-2), generalising the three-target / two-ligand docking analysis published there to a reusable many-to-many workflow. Hyun Ho Song (first author of this protocol) is a co-author of the published primary paper; the methodology presented here is therefore directly traceable to a peer-reviewed source.

The protocol's novelty rests on two engineering additions on top of vetted tools (ColabFold, P2Rank, GNINA 1.3, RDKit, PLIP, Snakemake): (i) an **automated pre-docking biological triage layer** (DeepTMHMM with optional Phobius + five rules) that removes integral-membrane and short-fragment ORFs before any GPU-bound work; and (ii) a **zero-parameter pocket-driven auto-box** that derives docking-box dimensions from the P2Rank pocket centre and the ligand maximum dimension. To our knowledge no published multi-target docking pipeline ships both. The full eight-tool comparison (DAD vs DiffDock, AutoDock Vina, GNINA standalone, AFPAP, AlphaFill, EquiBind, RFAA, Boltz-1) is given in Section 5 of the manuscript.

Current evidence is strong but deliberately bounded. The Tier-1 deterministic replay reproduces the supporting-paper docking values to numerical equality on six pairs (Δ = 0; 6 / 6 PASS). An open-license RCSB co-crystal seed of 16 cases redocked with GNINA 1.3.2 on a single RTX 2060 GPU achieves 16 / 16 GNINA execution, 12 / 16 (75 %) top-pose RMSD ≤ 2 Å, and 15 / 16 (94 %) best-of-9 search recovery; the CNN pose score classifies top-pose PASS / FAIL across cases with AUROC 0.958 but exhibits a within-case ranking limitation that we discuss explicitly in Limitations. A Tier-2 ~ 60-pair generalisation set across BindingDB / BioLiP2 / CrossDocked2020 / ChEMBL and a decoy / enrichment evaluation are **planned as Phase D / companion-study scope** and are not claimed in this manuscript.

DAD's reference-data dependency tree is entirely free and registration-free (RCSB CC0, BindingDB CC BY 3.0, BioLiP2, CrossDocked2020 CC0/MIT, ChEMBL CC BY-SA 4.0, AlphaFold DB CC BY 4.0); PDBbind and CASF-2016 are excluded by user policy. We believe this makes DAD an unusually clean fit for Nature Protocols' reproducibility expectation, since any reader can re-execute the entire validation pipeline without institutional gating.

The supporting primary paper is peer-reviewed and published in *iScience* 2026, satisfying Nature Protocols' supporting-research-paper requirement. The protocol manuscript is intentionally bounded to the current evidence: Tier 1 deterministic replay regression (six pairs from the published primary paper) and the 16-case open RCSB co-crystal seed external redocking. A live forward worked example from raw FASTA + SMILES through report outputs is queued as the next deliverable, and Tier-2 generalisation plus decoy / enrichment evaluation are signposted in the manuscript as Phase D / companion-study scope rather than claimed here.

Code is available under MIT (with the GNINA wrapper kept separate for upstream GPL-2.0 compliance) at https://github.com/hyunho-song09/dad-protocol (public, available at acceptance). The one-click Colab runner is at https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb. All quantitative claims are reproducible from `06_Report/Mr_Repro/results/benchmark/validation_table.tsv` and `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_results.tsv`. Five suggested reviewers are listed in `submission_metadata/suggested_reviewers.md`.

We thank the editors for their consideration.

Sincerely,
Do Yup Lee (corresponding author; <rome73@snu.ac.kr>)
On behalf of: Hyun Ho Song, Do Yup Lee.

---

*Word count target: ~ 400. Codex Round 6 unsafe-claim ban list is enforced (no "fully validated", "Tier 2 completed", "best-of-9 operational", or "exact identical scores"). Remaining placeholders for the camera-ready submission: Zenodo DOI (mint at acceptance of this protocol manuscript), submission date, exact iScience figure / supplementary identifier for the primary-paper Tier 1 ground truth (corresponding-author verification), and the iScience article id / DOI suffix from the published landing page; see `Publication/submission_metadata/PLACEHOLDER_RESOLUTION_GUIDE.md`.*
