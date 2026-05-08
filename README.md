# DAD: Dynamic Affinity Dock

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

DAD is a user-input protein-ligand docking notebook and code package.

## What It Does

- accepts user protein sequence or FASTA;
- accepts unnamed SMILES or `name:SMILES`;
- prepares 3D ligand SDF files;
- reuses cached/uploaded PDB files or predicts with automatic ESMFold/ColabFold AF2 fallback;
- creates an automatic docking box;
- runs GNINA docking when available;
- exports ranked results, plots, and a reproducibility JSON.

## Quick Start

1. Open the Colab notebook.
2. Run `0. Setup Environment`.
3. If condacolab restarts the kernel, run all cells again.
4. Paste your protein sequence in `custom_protein_fasta`.
5. Paste a SMILES string in `custom_ligand_smiles`.
6. Keep `STRUCTURE_MODE="auto"` to predict from sequence, or provide your own PDB with `existing_or_upload`.
7. Download `docking_results.tsv` or the result zip.


## Reuse Existing Structures

Keep the same `job_name` and change `custom_ligand_smiles` to score new ligands against the cached PDB. The notebook stores reusable structures in `WORK_DIR/structure_cache` or Drive `_structure_cache` when Drive output is enabled. Long proteins that exceed ESMFold API limits fall back to ColabFold AF2.

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
