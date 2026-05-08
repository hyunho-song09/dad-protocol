# DAD — Dynamic Affinity Dock

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**A reproducible Nature Protocols workflow for many-to-many protein–metabolite docking with pre-docking biological triage and zero-parameter defaults.**

DAD couples (i) automatic topology-aware ORF triage (DeepTMHMM with optional Phobius) and (ii) zero-parameter many-to-many GNINA scoring on a Snakemake DAG with three frozen runtimes (WSL2 / Colab / Docker–Apptainer). Tier 1 (six dipeptide × bacterial sensor cases) reproduces the supporting research paper *Sung et al. iScience 2026* with bit-equal accuracy; an independent 16-case open RCSB co-crystal seed achieves 12/16 top-pose RMSD ≤ 2.0 Å.

## Quick start (Colab, no install)

1. Click the **Open In Colab** badge above.
2. Set runtime to **GPU (T4)**.
3. **Run all cells**. Tier 1 replay completes in seconds; live RCSB redocking ≈ 4 min on T4.

Or clone locally:

```bash
git clone https://github.com/hyunho-song09/dad-protocol.git
cd dad-protocol/code
python -m pytest tests        # unit tests, no GPU
python ../source_code/fig04_cnn_roc.py   # regenerate Figure 4 from source TSV
```

## Two contributions

1. **Pre-docking biological triage** — five rules (length / signal-peptide clipping / TM topology / dock-region length / functional class) score each ORF before any structure or docking call. The triage output (`verdict ∈ {accept, downrank, exclude}` paired with biological status) gates the rest of the workflow and is, to our knowledge, the only such layer shipped natively in a published docking protocol (`Section 5, Table 2` of the manuscript).
2. **Zero-parameter many-to-many matrix** — N proteins × M ligands → an N×M GNINA score matrix, no parameters to tune. The protocol enforces `exhaustiveness = 32, num_modes = 9, seed = 0` and a deterministic auto-box rule (`side = max(22 Å, ligand_max_dim + 10 Å)`).

## What the package contains

| Folder | Purpose |
|--------|---------|
| `manuscript/` | DAD manuscript v1 (Markdown source for Nature Protocols submission) |
| `supplementary/` | Supplementary Information text + raw data TSVs |
| `figures/` | Five main + one SI figure as PDF / SVG / PNG dpi=300, with source code and source data |
| `tables/` | Six publication-ready tables (CSV) |
| `code/` | Cleaned Python package (`dad/`), Snakemake workflow, Colab/Jupyter notebooks, 134 unit tests |
| `references/` | `references.bib` (33 entries) + `citation_index.md` |
| `submission_metadata/` | Cover letter, suggested reviewers, data and code availability, conflict of interest, reporting summary, and the `USER_INPUTS_FILLED.yaml` master config |
| `source_data_audit.md` | Per-figure source TSV / generation script / SHA256 audit |
| `DAD_protocol.ipynb` | One-click Colab notebook (Colab Forms, hidden code, GUI inputs) |

## Citation

If you use DAD, please cite:

- **DAD protocol**: Song HH, Lee DY. DAD: Dynamic Affinity Dock — many-to-many protein–metabolite docking with pre-docking biological triage. *Nature Protocols*. 2026 (in preparation).
- **Supporting primary research paper**: Sung J-Y *et al.* Keratin degradation reflects a starvation survival strategy in *Fervidobacterium islandicum* AW-1. *iScience*. 2026; PII S2589-0042(26)01130-2. https://www.cell.com/iscience/fulltext/S2589-0042(26)01130-2

The full reference list and per-section citation index are in `references/`.

## License

DAD core (`code/dad/`) and the Snakemake workflow are released under **MIT**. Third-party dependencies (GNINA GPL-2.0, Open Babel GPL-2.0, PLIP GPL-2.0, ColabFold MIT, AlphaFold weights CC BY 4.0, RDKit BSD-3, Snakemake MIT, ChimeraX UCSF non-commercial, etc.) are listed with version pins in `submission_metadata/code_availability.md` and installed separately by the user. The DAD repository does **not** redistribute the GNINA binary.

Documentation files (`*.md`) are released under **CC BY 4.0**.

## Reproducibility

- **Frozen runtimes**: WSL2 + Ubuntu 24.04.3 (RTX 2060), Colab T4, Docker / Apptainer (HPC). Exact version pins: `code/workflow/runtime_freeze.md` (mirrored in `submission_metadata/code_availability.md`).
- **Source-data audit**: `source_data_audit.md` lists, per figure and table, the source TSV, the generation script, the generation date, and the SHA256 checksum.
- **Tests**: `code/tests/` runs without a GPU and verifies the Tier 1 replay regression and the 16-case redocking aggregation.

## Status

This repository accompanies a Nature Protocols Protocol manuscript currently in preparation (corresponding author: Do Yup Lee, Seoul National University, `rome73@snu.ac.kr`). The supporting primary research paper is *Sung et al. iScience 2026* (peer-reviewed published).
