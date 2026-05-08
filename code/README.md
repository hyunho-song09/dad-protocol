# DAD Code Package

[Open the Colab notebook](https://colab.research.google.com/github/hyunho-song09/dad-protocol/blob/main/DAD_protocol.ipynb)

## Two-Phase Workflow

| Phase | Steps | Input | Output |
|---|---|---|---|
| **Phase A** — Structure prep (run once) | §0–4 | Multi-FASTA | PDB + pocket cache, `structure_registry.tsv` |
| **Phase B** — Selective ligand scoring | §5–10 | Multi-SMILES + widget selection | `docking_master.csv`, heatmap |

## What It Does

| Stage | Description |
|---|---|
| Input (Phase A) | Multi-FASTA; raw sequences auto-named Protein_1, Protein_2, ... |
| Structure (Phase A) | `colabfold_batch` multi-FASTA CLI (default, `--num-models 1 --num-recycle 3 --sort-queries-by length`); or ingest AF3 Server CIF results; or ESMFold API; or user PDB |
| Pocket (Phase A) | P2Rank batch on all proteins in structure_registry |
| Input (Phase B) | Multi-SMILES (`name:SMILES` per line); SHA-keyed SDF cache |
| Selection (Phase B) | ipywidgets multi-select: choose protein subset × ligand subset |
| Docking (Phase B) | GNINA cross-product of selected pairs only; cache-hit pairs skipped |
| Export | `docking_master.csv` (append-only), ranked table, CNN-affinity heatmap |

## User Workflow

**Phase A (run once per protein set):**
1. §1 Paste multi-FASTA sequences.
2. §2 Set `STRUCTURE_MODE = "colabfold_af2"` (default) or `"af3_results"` if you have AlphaFold Server results.
3. §3–4 P2Rank pocket detection + `structure_registry.tsv` summary.

**Phase B (repeat with new SMILES or new subset):**
4. §5 Paste SMILES (`name:SMILES` per line).
5. §6 Select proteins and ligands in the widget UI.
6. §7 Click **Run docking on selection**.
7. §8 Export `docking_master.csv`.

Re-running Phase B with new SMILES does **not** re-run Phase A.

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
