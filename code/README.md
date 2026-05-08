# DAD: Dynamic Affinity Dock

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/[USER_TO_PROVIDE_GITHUB_OWNER]/[USER_TO_PROVIDE_GITHUB_REPO]/blob/main/code/notebooks/DAD_protocol.ipynb)

**Google Colab에서 바로 실행 / Run directly in Google Colab:**
위 배지를 클릭하면 브라우저에서 DAD 프로토콜이 즉시 열립니다. GPU 런타임(T4 권장)을 선택한 후 "모두 실행"을 누르면 Tier 1 재현 분석(3 proteins × 2 ligands)이 GPU 없이 약 2분 안에 완료됩니다.

Click the badge above to open the DAD protocol notebook in your browser. Select a GPU runtime (T4 recommended), then click "Run all". The default Tier 1 replay mode completes in under 2 minutes without GPU. For live GNINA docking, switch `exec_mode` to `"live"` in §1 Input.

---

DAD is a Nature Protocols-style many-to-many protein-metabolite docking workflow for bacterial chemosensory receptors and substrate-binding proteins.

The protocol contribution is not just another docking wrapper. It combines:

- pre-docking biological triage of ORFs using membrane/topology evidence and rules R1-R5;
- zero-parameter pocket box generation for many protein-ligand pairs;
- replay, live redocking, figure source-data, and manifest outputs for reproducibility.

## Current Validation Status

RCSB seed benchmark: 16 protein-ligand co-crystal cases.

- Top-pose success: 12/16 cases with RMSD <= 2.0 A.
- Best-of-9 success: 15/16 cases with RMSD <= 2.0 A.
- CNN score classifier: AUROC 0.958 for top-pose PASS/FAIL.
- Tier 1 replay benchmark: expected 6/6 PASS.

These numbers describe the current publication package draft, not a final journal claim until the manuscript tables, captions, and source data are locked together.

## Directory Structure

```text
Publication/code/
  dad/                    Python package source
  workflow/               Snakemake workflow and zero-parameter config
  notebooks/              Colab/user notebooks
  tests/                  Pytest regression tests
  README.md               This file
  CITATION.cff            Citation metadata
```

## Quick Start: No-GPU Replay Smoke Test

From the project root:

```powershell
cd d:\project\experiment\DAD
py -3 06_Report\Mr_Repro\benchmark\run_replay.py
```

Expected result:

- `06_Report/Mr_Repro/replay_results/validation_table.tsv`
- `06_Report/Mr_Repro/replay_results/docking_master.csv`
- `06_Report/Mr_Repro/replay_results/manifest.json`
- summary: 6/6 Tier 1 cases PASS

## Quick Start: Import the Publication Code

Install the publication snapshot in editable mode:

```powershell
cd d:\project\experiment\DAD
py -3 -m pip install -e Publication\code --no-deps
py -3 -c "import dad; print(dad.__version__)"
```

If you do not want to modify the active Python environment, set `PYTHONPATH` instead:

```powershell
cd d:\project\experiment\DAD
$env:PYTHONPATH = "$PWD\Publication\code"
py -3 -c "import dad; print(dad.__version__)"
```

## Snakemake Dry Run

```powershell
cd d:\project\experiment\DAD\Publication\code
$env:PYTHONPATH = "$PWD"
snakemake --snakefile workflow\Snakefile --configfile workflow\config.yaml --dryrun --printshellcmds
```

The workflow resolves paths back to the project root when it is run from the publication snapshot inside this repository. For an archived release, preserve the root-level `01_Tier1_input` and `06_Report` directories next to `Publication`.

## Pytest Regression Tests

Recommended environment:

```bash
conda create -n dad-lite -c conda-forge python=3.11 rdkit biopython numpy scipy pandas pytest
conda activate dad-lite
cd /path/to/DAD/Publication/code
pytest tests -q
```

Windows PowerShell path:

```powershell
cd d:\project\experiment\DAD\Publication\code
pytest tests -q --basetemp=..\pytest_tmp
```

## Colab Entry Point

Open either notebook:

- `DAD_protocol.ipynb` at the repository root
- `Publication/code/notebooks/DAD_protocol.ipynb`

Use a T4 GPU runtime for live GNINA/ColabFold execution. CPU is enough for replay-only checks and static source-data regeneration.

## Dependencies

Runtime Python dependencies: Python >= 3.10, numpy, scipy, scikit-learn, matplotlib, rdkit, biopython, pandas.

External tools downloaded or installed separately:

| Tool | Version used in package | License | Note |
|---|---:|---|---|
| GNINA | 1.3.2 | GPL-2.0 | Not redistributed; downloaded at runtime |
| ColabFold | 1.6.1 | Apache-2.0 | Used through Colab/local installation |
| P2Rank | 2.4.2 | MIT | External executable |
| DeepTMHMM | latest service/package | MIT | Used for topology triage |

## License

The `dad/` Python package is released under the MIT License. Third-party tools and datasets retain their own licenses; see the manuscript Data Availability and Code Availability files.
