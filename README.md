# DAD: Dynamic Affinity Dock

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

DAD is a user-input protein-ligand docking notebook and code package.

## What It Does

- accepts multi-FASTA protein input;
- accepts unnamed SMILES, `name:SMILES`, or semicolon-separated SMILES;
- prepares 3D ligand SDF files;
- separates Phase A structure preparation from Phase B ligand scoring;
- reuses cached structures and selected protein-ligand pair outputs;
- runs GNINA docking when available.

## Quick Start

1. Open the Colab notebook.
2. Run `0. Setup Environment`.
3. If condacolab restarts the kernel, run all cells again.
4. In SS1, paste multi-FASTA sequences.
5. In SS2, keep `STRUCTURE_MODE="colabfold_af2"` or choose `auto`, `af3_results`, `esmfold_api`, or `user_pdb`.
6. In SS5, paste SMILES entries in `smiles_text`.
7. Select protein-ligand pairs and export `phase_b/docking_master.csv`.


## Two-Phase Reuse

Phase A writes `phase_a/structure_registry.tsv` and reusable PDB files in `phase_a/structure_cache`. Phase B can be rerun with new `smiles_text` without re-running structure prediction.

## Input Examples

Raw protein sequence:

```text
MRNMSIFMKVMVIVLILALGMIVIGVYSTFAL...
```

FASTA:

```text
>ProteinA
MRNMSIFMKVMVIVLILALGMIVIGVYSTFAL...
```

Unnamed ligand:

```text
CC[C@H](C)[C@@H](C(=O)O)NC(=O)[C@H](C)N
```

Named ligand:

```text
Ala-Ile:CC[C@H](C)[C@@H](C(=O)O)NC(=O)[C@H](C)N
```

## Repository Map

| Path | Purpose |
|---|---|
| `DAD_protocol.ipynb` | main Colab notebook |
| `code/dad/` | Python package |
| `code/workflow/` | Snakemake workflow |
| `code/tests/` | unit tests |
| `figures/`, `tables/`, `supplementary/` | reference examples and documentation assets |

## Local Package Check

```bash
git clone https://github.com/hyunho-song09/dad-protocol.git
cd dad-protocol/code
python -m pip install -e . --no-deps
python -c "import dad; print(dad.__version__)"
```

## License

DAD core code is MIT licensed. Documentation is CC BY 4.0. External tools keep their own licenses. This repository does not redistribute the GNINA binary.
