# DAD Colab Usage Guide

## Quick Start

1. Open `DAD_protocol.ipynb`.
2. Run `0. Setup Environment`.
3. If the kernel restarts after condacolab, run all cells again.
4. Paste your protein sequence or FASTA.
5. Paste a ligand SMILES string.
6. Keep `STRUCTURE_MODE="esmfold_api"` to predict from sequence, or switch to `existing_or_upload` for your own PDB.
7. Run docking and download the result archive.

## Progress Bars

Every notebook step prints a progress bar:

```text
Setup Environment: [########------------] 50% 4/8 chemistry packages
```

Long-running commands such as ESMFold API prediction and GNINA still report only step-level progress.

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
| `esmfold_api` | Default sequence-to-PDB prediction without AlphaFold/TensorFlow installs |
| `existing_or_upload` | You already have PDB files |
| `colabfold` | Optional advanced mode on Python versions compatible with AlphaFold/TensorFlow |

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

- keep `STRUCTURE_MODE="esmfold_api"` and retry;
- upload `Protein_1.pdb` to `WORK_DIR/structures`;
- set `EXISTING_PDB_DIR`.

### GNINA not found

Run Setup with:

```python
INSTALL_GNINA = True
```

### GPU quota exceeded

Use `exec_mode="prepare_only"` to prepare inputs without docking.
