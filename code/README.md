# DAD Code Package

[Open the Colab notebook](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)

## User Workflow

1. Paste protein sequence or FASTA.
2. Paste SMILES or `name:SMILES`.
3. Provide a PDB or run ColabFold.
4. Generate ligand SDF files.
5. Run GNINA.
6. Export ranked TSV results.

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
