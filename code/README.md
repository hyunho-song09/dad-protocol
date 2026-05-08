# DAD Code Package

[Open the Colab notebook](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)

## Two-Phase Workflow

| Phase | Steps | Input | Output |
|---|---|---|---|
| Phase A: structure prep | SS0-SS4 | Multi-FASTA | PDB cache, pocket cache, `structure_registry.tsv` |
| Phase B: selective ligand scoring | SS5-SS10 | Multi-SMILES and selected pairs | `docking_master.csv`, ranked table, heatmap |

## What It Does

| Stage | Description |
|---|---|
| Input | Multi-FASTA; raw sequences are auto-named `Protein_1`, `Protein_2`, ... |
| Structure | Default `colabfold_af2`; alternatives are `af3_results`, `esmfold_api`, `user_pdb`, and `auto` |
| Pocket | P2Rank batch when available; protein-center fallback when P2Rank is unavailable |
| Ligand | Multi-SMILES, one `name:SMILES` per line or semicolon-separated; SHA-keyed SDF cache |
| Selection | `ipywidgets` multi-select for protein subset and ligand subset |
| Docking | GNINA cross-product of selected pairs only; completed pair cache is reused |
| Export | `phase_b/docking_master.csv`, ranked table, and diagnostic heatmap |

## User Workflow

Phase A, run once per protein set:

1. SS1: paste multi-FASTA sequences.
2. SS2: keep `STRUCTURE_MODE="colabfold_af2"` or choose another backend.
3. SS3-SS4: run pocket detection and write `structure_registry.tsv`.

Phase B, repeat as needed:

1. SS5: paste new SMILES.
2. SS6: select proteins and ligands.
3. SS7: run docking for selected pairs.
4. SS8: update `docking_master.csv`.

Re-running Phase B with new SMILES does not re-run Phase A.

## Install

```bash
cd dad-protocol/code
python -m pip install -e . --no-deps
python -c "import dad; print(dad.__version__)"
```

## Tests

```bash
conda create -n dad-lite -c conda-forge python=3.11 rdkit biopython numpy scipy pandas pytest
conda activate dad-lite
pytest tests -q
```

## Key Files

| Path | Purpose |
|---|---|
| `dad/` | Python package |
| `workflow/` | Snakemake workflow |
| `notebooks/DAD_protocol.ipynb` | notebook copy |
| `tests/` | pytest suite |
| `pyproject.toml` | editable install metadata |

## External Tools

GNINA, ColabFold, P2Rank, RDKit, and Biopython keep their own licenses and installation rules. The repository does not redistribute GNINA.
