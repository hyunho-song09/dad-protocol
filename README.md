# DAD: Dynamic Affinity Dock

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

DAD is a reproducible many-to-many protein-metabolite docking workflow for bacterial receptor studies.

Core features:

- topology-aware ORF triage before docking;
- zero-parameter pocket-based docking boxes;
- Tier 1 replay mode for no-GPU testing;
- live GNINA docking mode;
- source-data, figures, tables, and reproducibility metadata.

## Quick Start

### Colab

1. Click `Open In Colab`.
2. Click `Runtime > Run all`.
3. If condacolab restarts the kernel, click `Runtime > Run all` again.
4. Keep `exec_mode = "tier1_replay"` for the fastest smoke test.

### Local Tests

```bash
git clone https://github.com/hyunho-song09/dad-protocol.git
cd dad-protocol/code
conda create -n dad-lite -c conda-forge python=3.11 rdkit biopython numpy scipy pandas pytest
conda activate dad-lite
pytest tests -q
```

## Validation Snapshot

- Tier 1 replay: 6/6 expected PASS.
- RCSB seed: 16 cases.
- Top-pose success: 12/16 at RMSD <= 2.0 A.
- Best-of-9 success: 15/16 at RMSD <= 2.0 A.
- CNN score AUROC: 0.958 for top-pose PASS/FAIL.

## Repository Map

| Path | Purpose |
|---|---|
| `DAD_protocol.ipynb` | one-click Colab notebook |
| `code/dad/` | Python package |
| `code/workflow/` | Snakemake workflow |
| `code/tests/` | unit tests |
| `figures/` | figures, source code, source data |
| `tables/` | publication tables |
| `supplementary/` | supplementary text and data |
| `submission_metadata/` | submission support files |

## License

DAD core code is MIT licensed. Documentation is CC BY 4.0. External tools keep their own licenses. This repository does not redistribute the GNINA binary.
