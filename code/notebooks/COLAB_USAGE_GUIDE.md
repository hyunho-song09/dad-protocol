# DAD Colab Usage Guide

## Quick Start

1. Open `DAD_protocol.ipynb`.
2. Run `0. Setup Environment`.
3. If the kernel restarts after condacolab, run all cells again.
4. Paste your protein sequence or FASTA.
5. Paste a ligand SMILES string.
6. Keep `STRUCTURE_MODE="colabfold"` to predict from sequence, or switch to `existing_or_upload` for your own PDB.
7. Run docking and download the result archive.

## Progress Bars

Every notebook step prints a progress bar:

```text
Setup Environment: [########------------] 50% 4/8 chemistry packages
```

Long-running commands such as ColabFold and GNINA still report only step-level progress.

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
| `existing_or_upload` | You already have PDB/CIF files |
| `colabfold` | You want ColabFold to predict structures from FASTA |

For `existing_or_upload`, put files in:

```text
WORK_DIR/structures/
```

or set:

```text
EXISTING_PDB_DIR=/path/to/pdbs
```

The expected file name is `Protein_1.pdb` unless the FASTA header gives another ID.

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

- upload `Protein_1.pdb` to `WORK_DIR/structures`;
- set `EXISTING_PDB_DIR`;
- set `STRUCTURE_MODE="colabfold"`.

### GNINA not found

Run Setup with:

```python
INSTALL_GNINA = True
```

### GPU quota exceeded

Use `exec_mode="prepare_only"` to prepare inputs without docking.
