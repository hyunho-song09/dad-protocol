# DAD Workflow: Snakemake Usage Guide

**Files:** `workflow/Snakefile`, `workflow/config.yaml`  
**Snapshot:** publication package, 2026-05-08

## Prerequisites

```bash
pip install "snakemake>=7.0"
```

For local execution, make the source package importable:

```powershell
cd d:\project\experiment\DAD\Publication\code
py -3 -m pip install -e . --no-deps
```

No-install alternative:

```powershell
cd d:\project\experiment\DAD\Publication\code
$env:PYTHONPATH = "$PWD"
```

For conda/Linux:

```bash
cd /path/to/DAD/Publication/code
export PYTHONPATH="$PWD"
```

Optional live-mode tools:

- GNINA v1.3.2
- P2Rank v2.4.2
- DeepTMHMM or Phobius-compatible topology input
- ColabFold/AlphaFold structure backend

## Dry Run

From the publication code directory:

```powershell
cd d:\project\experiment\DAD\Publication\code
$env:PYTHONPATH = "$PWD"
snakemake --snakefile workflow\Snakefile --configfile workflow\config.yaml --dryrun --printshellcmds
```

The Snakefile resolves config data paths back to the project root by locating the root-level `06_Report` and `01_Tier1_input` directories. This lets Claude or a reviewer run the publication snapshot without copying the workflow back into `06_Report/Mr_Pipeline`.

## Execution Modes

The `execution.mode` key in `workflow/config.yaml` controls the branch.

| Mode | Description | GPU required |
|---|---|---|
| `tier1_replay` | Uses AW1_ref/precomputed assets and ground-truth scores. Validates data flow without GNINA or ColabFold. | No |
| `live` | Runs structure/pocket/docking tools where configured. | Yes for practical full runs |

Override on the command line:

```powershell
snakemake --snakefile workflow\Snakefile --configfile workflow\config.yaml --cores 4 --config execution.mode=tier1_replay
```

## Full Tier 1 Replay

Prefer the dedicated replay harness for reviewer smoke tests:

```powershell
cd d:\project\experiment\DAD
py -3 06_Report\Mr_Repro\benchmark\run_replay.py
```

Expected outputs:

- `06_Report/Mr_Repro/replay_results/validation_table.tsv`
- `06_Report/Mr_Repro/replay_results/docking_master.csv`
- `06_Report/Mr_Repro/replay_results/manifest.json`
- expected result: 6/6 PASS

## Live Redocking

Live redocking is better run from WSL/Linux or Colab because GNINA is GPU-oriented.

```bash
python 06_Report/Mr_Repro/benchmark/run_rcsb_redocking.py \
  --gnina 06_Report/Mr_Repro/tools/gnina/run_gnina_wsl.sh \
  --exhaustiveness 32 \
  --num-modes 9 \
  --seed 0
```

For a quick smoke test, add one of:

```bash
--limit 2
--case-id CRP_CMP_1I5Z
```

Current seed benchmark expectation:

- 16/16 GNINA jobs complete
- 12/16 top poses PASS at RMSD <= 2.0 A
- 15/16 best-of-9 poses PASS at RMSD <= 2.0 A

## Troubleshooting

- If `ModuleNotFoundError: dad` appears, set `PYTHONPATH` to `Publication/code`.
- If GNINA fails on Windows, use WSL/Linux or Colab with GPU runtime.
- If a full Snakemake live run is too slow, first run `run_replay.py`, then run `run_rcsb_redocking.py --limit 2` as a smoke test.
