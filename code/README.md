# DAD Code Package

[Open the Colab notebook](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)

## What It Does

| Stage | Description |
|---|---|
| Input | Accepts raw protein sequences or FASTA; raw SMILES or `name:SMILES` |
| Triage | Filters by length, signal-peptide clipping, TM-helix count, dock-region length, and functional class |
| Structure | Accepts AlphaFold 3 results (default; from alphafoldserver.com or local install) or, optionally, predicts via ESMFold API or ColabFold AF2 |
| Pocket | P2Rank binding-site prediction |
| Docking | GNINA scoring with per-pose CNN affinity |
| Export | Ranked TSV with ipTM / pTM / pLDDT / CNN-affinity columns |

## User Workflow

1. Paste protein sequence or FASTA.
2. Paste SMILES or `name:SMILES`.
3. Run triage.
4. Triage done — run structure step.
5a. (Preferred) Download AF3 results from alphafoldserver.com, set `STRUCTURE_MODE = "af3_results"`, enter the folder path in `AF3_RESULTS_PATH`.
5b. (Alternative) Set `STRUCTURE_MODE` to `esmfold_api` or `colabfold_af2` for in-notebook prediction.
6. Run GNINA.
7. Export ranked TSV results.

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
