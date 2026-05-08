# DAD Code Package

[Open the Colab notebook](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)

DAD is a protein-metabolite docking workflow with:

- topology-aware ORF triage;
- zero-parameter pocket box generation;
- Tier 1 replay mode;
- live GNINA docking mode;
- figure/source-data regeneration.

## Quick Start

```bash
git clone https://github.com/hyunho-song09/dad-protocol.git
cd dad-protocol/code
python -m pip install -e . --no-deps
python -c "import dad; print(dad.__version__)"
```

## Colab

1. Open `DAD_protocol.ipynb`.
2. Run all cells.
3. If Setup installs condacolab, the kernel restarts.
4. After restart, run all cells again.
5. Use `tier1_replay` for the fast no-GPU test.
6. Use `live` only when GNINA/ColabFold execution is needed.

## Tests

```bash
conda create -n dad-lite -c conda-forge python=3.11 rdkit biopython numpy scipy pandas pytest
conda activate dad-lite
cd dad-protocol/code
pytest tests -q
```

## Workflow Dry Run

```bash
cd dad-protocol/code
export PYTHONPATH="$PWD"
snakemake --snakefile workflow/Snakefile --configfile workflow/config.yaml --dryrun --printshellcmds
```

## Validation Snapshot

- RCSB seed: 16 cases.
- Top-pose success: 12/16 at RMSD <= 2.0 A.
- Best-of-9 success: 15/16 at RMSD <= 2.0 A.
- CNN score AUROC: 0.958 for top-pose PASS/FAIL.
- Tier 1 replay: expected 6/6 PASS.

## Key Files

| Path | Purpose |
|---|---|
| `dad/` | Python package |
| `workflow/` | Snakemake workflow |
| `notebooks/DAD_protocol.ipynb` | Colab notebook copy |
| `tests/` | pytest suite |
| `pyproject.toml` | editable install metadata |

## External Tools

GNINA, ColabFold, P2Rank, DeepTMHMM/Phobius, RDKit, and Biopython keep their own licenses and installation rules. The repository does not redistribute the GNINA binary.
