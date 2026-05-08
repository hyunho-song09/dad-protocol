# DAD Colab Usage Guide

## Quick Start

1. Open `DAD_protocol.ipynb`.
2. Run `0. Setup Environment`.
3. If condacolab restarts the kernel, run all cells again.
4. In SS1, paste one or more protein FASTA records.
5. In SS2, keep `STRUCTURE_MODE="colabfold_af2"` for direct AF2 prediction, or choose `auto`, `af3_results`, `esmfold_api`, or `user_pdb`.
6. In SS3-SS4, create `phase_a/structure_registry.tsv`.
7. In SS5, paste one or more SMILES entries.
8. In SS6-SS9, select protein-ligand pairs, run GNINA, aggregate results, and export pose visualizations.

## Structure Modes

| Option | Use When |
|---|---|
| `colabfold_af2` | Default ColabFold AF2 prediction with multi-FASTA batching |
| `auto` | Try cache, then ESMFold API, then ColabFold AF2 fallback |
| `af3_results` | Ingest pre-computed AlphaFold Server or local AF3 CIF results |
| `esmfold_api` | ESMFold API for short proteins, with optional ColabFold fallback |
| `user_pdb` | Use uploaded PDBs or files in `USER_PDB_DIR` |

## Reuse Mode

Phase A is sequence-only and ligand-independent. To test new ligands against the same proteins:

1. Keep the same work directory.
2. Skip SS1-SS4 if `phase_a/structure_registry.tsv` already exists.
3. Run SS5 onward with new `smiles_text`.

The notebook reloads `structure_registry.tsv` when Phase B starts in a fresh runtime.

## Accepted Inputs

FASTA:

```text
>ProteinA
MRNMSIFMKVMVIVLILALGMIVIGVYSTFAL...
>ProteinB
MNNNKQQQQ...
```

SMILES:

```text
Ala-Ile:CC[C@H](C)[C@@H](C(=O)O)NC(=O)[C@H](C)N
LigandB:CCO
```

Semicolon-separated SMILES are also accepted:

```text
LigandA:CCO;LigandB:CCN
```

## Outputs

| File | Content |
|---|---|
| `phase_a/structure_registry.tsv` | protein structure paths, confidence, pocket metadata |
| `phase_a/structure_cache/*.pdb` | reusable PDB files keyed by sequence SHA |
| `phase_a/pocket_cache/` | P2Rank outputs when available |
| `phase_b/ligands/*.sdf` | reusable 3D ligand files keyed by SMILES SHA |
| `phase_b/runs/<case_id>/docked.sdf` | GNINA output per selected pair |
| `phase_b/docking_master.csv` | append-only selected-pair result table |
| `phase_b/visualizations/*_contacts.tsv` | contact residues within the selected distance cutoff |
| `phase_b/visualizations/*_complex.pdb` | merged receptor-ligand complex for external viewers |
| `phase_b/visualizations/*.cxc` | ChimeraX styling script with contact labels and PNG export command |
| `phase_b/visualizations/interaction_contacts.tsv` | combined contact-residue table |
| `manifest.json` | reproducibility footprint |

## Troubleshooting

### ESMFold HTTP 413

HTTP 413 means the ESMFold API rejected the sequence size. Use `STRUCTURE_MODE="auto"` or `STRUCTURE_MODE="colabfold_af2"`.

### ColabFold TensorFlow Import Error

The notebook applies the TensorFlow shared-object cleanup used in the reference AlphaFold2 notebook before importing ColabFold. Check:

```text
phase_a/tensorflow_colabfold_patch.log
phase_a/colabfold_install.log
phase_a/colabfold.log
```

### P2Rank Not Found

Pocket detection falls back to the protein geometric center, so GNINA can still run. Install P2Rank for better pocket centers.

### GNINA Not Found

Run SS0 again in a Linux or Colab runtime. If GNINA is still unavailable, SS7 reports `gnina_missing` for selected pairs.
