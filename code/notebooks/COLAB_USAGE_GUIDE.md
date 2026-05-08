# DAD Colab Usage Guide

## Quick Start

1. Open `DAD_protocol.ipynb`.
2. Run `0. Setup Environment`.
3. If the kernel restarts after condacolab, run all cells again.
4. Paste your protein sequence or FASTA.
5. Paste a ligand SMILES string.
6. Keep `STRUCTURE_MODE="auto"` to predict from sequence, or switch to `existing_or_upload` for your own PDB.
7. Run docking and download the result archive.

## Progress Bars

Every notebook step prints a progress bar:

```text
Setup Environment: [########------------] 50% 4/8 chemistry packages
```

Long-running commands such as ESMFold API, ColabFold AF2, and GNINA still report only step-level progress.

## Accepted Inputs

Protein sequence:

```text
MRNMSIFMKVMVIVLILALGMIVIGVYSTFAL...
```

FASTA:

```text
>ProteinA
MRNMSIFMKVMVIVLILALGMIVIGVYSTFAL...
```

SMILES:

```text
CC[C@H](C)[C@@H](C(=O)O)NC(=O)[C@H](C)N
```

Named SMILES:

```text
Ala-Ile:CC[C@H](C)[C@@H](C(=O)O)NC(=O)[C@H](C)N
```

Multiple ligands:

```text
LigandA:CCO;LigandB:CCN
```

## Structure Options

| Option | Use When |
|---|---|
| `auto` | Default: cache/uploaded PDB, then ESMFold API, then ColabFold AF2 fallback |
| `existing_or_upload` | You already have PDB files |
| `esmfold_api` | ESMFold API only, with optional ColabFold fallback |
| `colabfold` | Direct ColabFold AF2 using the AlphaFold2.ipynb backend |

For `existing_or_upload`, put files in:

```text
WORK_DIR/structures/
```

or set:

```text
EXISTING_PDB_DIR=/path/to/pdbs
```

The expected file name is `Protein_1.pdb` unless the FASTA header gives another ID.


## Reuse Mode

To test new ligands against the same protein:

1. Keep the same `job_name`.
2. Change only `custom_ligand_smiles`.
3. Run from `1. Input Data` onward.

The structure cell reuses `structure_manifest.tsv` and the sequence-hash cache before any new prediction.

## Outputs

| File | Content |
|---|---|
| `results/docking_results.tsv` | ranked docking results |
| `results/tables/ranked_results.tsv` | same ranked table |
| `results/tables/protein_summary.tsv` | best score per protein |
| `results/tables/ligand_summary.tsv` | best score per ligand |
| `results/figures/*.png` | result plots |
| `results/reproducibility_footprint.json` | runtime metadata |
| `<job_name>_results.zip` | result archive |

## Troubleshooting

### No ligands parsed

Paste either raw SMILES or `name:SMILES`. Raw SMILES is now accepted.

### No PDB structures found

Use one of:

- keep `STRUCTURE_MODE="auto"` and retry;
- upload `Protein_1.pdb` to `WORK_DIR/structures`;
- set `EXISTING_PDB_DIR`.


### ESMFold HTTP 413

HTTP 413 means the ESMFold API rejected the sequence size. Keep `STRUCTURE_MODE="auto"`; the notebook will switch to ColabFold AF2 and apply the TensorFlow crash patch used in the AlphaFold2.ipynb reference notebook.

For the most conservative AF2 setting, use `AF2_PRESET="high_accuracy"`. If you also set `AF2_RELAX_TOP_N=1`, the notebook installs OpenMM/PDBFixer before AMBER relaxation.

### GNINA not found

Run Setup with:

```python
INSTALL_GNINA = True
```

### GPU quota exceeded

Use `exec_mode="prepare_only"` to prepare inputs without docking.
