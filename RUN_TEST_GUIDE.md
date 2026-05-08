# DAD User Run and Test Guide

This guide answers where a reviewer or collaborator can run the DAD platform and what a successful smoke test should look like.

## 1. Colab Entry Points

Use one of these notebooks:

- `DAD_protocol.ipynb` at the repository root
- `Publication/code/notebooks/DAD_protocol.ipynb`

Recommended Colab setup:

1. Upload or mount the complete `DAD` folder.
2. Open the notebook.
3. Use `Runtime > Change runtime type > T4 GPU` for live GNINA/ColabFold execution.
4. Run setup cells first, then choose replay or live mode.

CPU Colab is acceptable only for replay/static checks. Use GPU for live redocking.

## 2. Fast Local Replay Test: No GPU

From Windows PowerShell:

```powershell
cd d:\project\experiment\DAD
py -3 06_Report\Mr_Repro\benchmark\run_replay.py
```

Expected outputs:

- `06_Report/Mr_Repro/replay_results/validation_table.tsv`
- `06_Report/Mr_Repro/replay_results/docking_master.csv`
- `06_Report/Mr_Repro/replay_results/manifest.json`

Expected result:

- 8 AW1 assets loaded
- 6 Tier 1 ground-truth cases loaded
- 6/6 validation cases PASS

This is the best first test for Claude or a reviewer because it does not require GPU, GNINA, ColabFold, or network access.

## 3. Publication Code Import Test

Editable install:

```powershell
cd d:\project\experiment\DAD
py -3 -m pip install -e Publication\code --no-deps
py -3 -c "import dad; print(dad.__version__)"
```

No-install alternative:

```powershell
cd d:\project\experiment\DAD
$env:PYTHONPATH = "$PWD\Publication\code"
py -3 -c "import dad; print(dad.__version__)"
```

## 4. Pytest Regression Tests

Recommended conda environment:

```bash
conda create -n dad-lite -c conda-forge python=3.11 rdkit biopython numpy scipy pandas pytest
conda activate dad-lite
cd /path/to/DAD/Publication/code
pytest tests -q
```

Windows path:

```powershell
cd d:\project\experiment\DAD\Publication\code
pytest tests -q --basetemp=..\pytest_tmp
```

Notes:

- Tests mock external tool calls.
- No GPU is required.
- RDKit is easiest through conda-forge.

## 5. Snakemake Dry Run

```powershell
cd d:\project\experiment\DAD\Publication\code
$env:PYTHONPATH = "$PWD"
snakemake --snakefile workflow\Snakefile --configfile workflow\config.yaml --dryrun --printshellcmds
```

The publication Snakefile now resolves paths back to the project root by looking for `06_Report` and `01_Tier1_input`.

## 6. Live RCSB Redocking

Use WSL/Linux or Colab GPU. Full Windows-native GNINA execution is not the recommended path.

WSL example:

```powershell
wsl -d Ubuntu -- bash -lc "cd /mnt/d/project/experiment/DAD && python3 06_Report/Mr_Repro/benchmark/run_rcsb_redocking.py --gnina 06_Report/Mr_Repro/tools/gnina/run_gnina_wsl.sh --exhaustiveness 32 --num-modes 9 --seed 0"
```

Quick smoke variants:

```powershell
wsl -d Ubuntu -- bash -lc "cd /mnt/d/project/experiment/DAD && python3 06_Report/Mr_Repro/benchmark/run_rcsb_redocking.py --gnina 06_Report/Mr_Repro/tools/gnina/run_gnina_wsl.sh --limit 2"
```

```powershell
wsl -d Ubuntu -- bash -lc "cd /mnt/d/project/experiment/DAD && python3 06_Report/Mr_Repro/benchmark/run_rcsb_redocking.py --gnina 06_Report/Mr_Repro/tools/gnina/run_gnina_wsl.sh --case-id CRP_CMP_1I5Z"
```

Expected current seed benchmark:

- 16/16 redocking jobs complete
- 12/16 top-pose RMSD <= 2.0 A
- 15/16 best-of-9 RMSD <= 2.0 A
- median top-pose RMSD about 0.56 A

## 7. Figure and Source-Data Regeneration

```powershell
cd d:\project\experiment\DAD
py -3 Publication\figures\source_code\fig04_cnn_roc.py
py -3 Publication\figures\source_data\generate_source_data.py
```

Expected files:

- `Publication/figures/source_data/fig04_source_data.csv`
- `Publication/figures/source_data/fig04_metrics.csv`
- `Publication/figures/source_data/fig03_source_data.csv`
- `Publication/figures/source_data/figS01_source_data.csv`
- `Publication/figures/source_data/README.md`
