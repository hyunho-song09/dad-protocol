# Code Availability Statement

> Mr_Prof, 2026-05-08. Drop-in for Nature Protocols Reporting Summary §"Code Availability".

---

All code underlying this protocol is open-source and available through the DAD project repository. The protocol's licensing policy keeps the GNINA wrapper separate from the DAD core in order to comply with upstream GPL-2.0 obligations.

## Repository

- **URL**: https://github.com/hyunho-song09/dad-protocol (corresponding author has completed the GitHub OAuth step; team-lead has proposed `dad-protocol` as a public repository — final URL pending)
- **Archived release**: Zenodo DOI on hold per corresponding author (2026-05-08); to be minted at acceptance of this protocol manuscript.
- **Git tag at submission**: `v0.1.0` (proposed)

## Licensing

| Component | License | Reason |
|-----------|---------|--------|
| `dad/core/*.py` (DAD library) | **MIT** | Permissive; allows commercial reuse |
| `06_Report/Mr_Pipeline/Snakefile` and rules | **MIT** | Same |
| `06_Report/Mr_Repro/benchmark/*.py` (replay + RCSB redocking runners) | **MIT** | Same |
| `06_Report/Mr_Repro/tools/gnina/run_gnina_wsl.sh` (GNINA wrapper) | **isolated**; wrapper is MIT, **but invokes the upstream GNINA binary which is GPL-2.0** | The wrapper is a thin CLI shim that sets `LD_LIBRARY_PATH` and execs the GNINA binary via `exec`. It does not statically link against GNINA. Users who redistribute the bundled GNINA binary must comply with GNINA's GPL-2.0 terms separately; the DAD repository itself does **not** redistribute the GNINA binary — `tools/gnina/run_gnina_wsl.sh` resolves the binary at run time from a user-installed path. |
| `06_Report/Mr_Repro/container-hpc/Dockerfile` | **MIT** for the Dockerfile itself; bundled GNINA inside the resulting image is GPL-2.0 | Users building the image should mark its derived form GPL-2.0 in their own redistribution. |
| Notebooks (`DAD_protocol.ipynb`, `Mr_Pipeline/colab/dad_run.ipynb`) | **MIT** | Same |
| Documentation (`*.md`) | **CC BY 4.0** | Standard scientific-text licence |

## Third-party software dependencies

The DAD pipeline depends on the following external tools, each under its own license. Users install these separately; DAD does not redistribute them.

| Tool | Version pin | License |
|------|-------------|---------|
| ColabFold (incl. AlphaFold 2) | ≥ 1.6.1 | MIT (ColabFold) + Apache 2.0 (AlphaFold) |
| GNINA | 1.3.2 (`master:f23dd2b`) | **GPL-2.0** |
| AutoDock Vina (via GNINA) | bundled in GNINA | Apache 2.0 |
| P2Rank / PrankWeb | ≥ 2.4 | MIT |
| DeepTMHMM | v1.0.24 | MIT |
| Phobius | (opt-in) | Academic free; commercial restricted |
| RDKit | 2024.09+ | BSD-3-Clause |
| Open Babel | 3.1+ | GPL-2.0 |
| PLIP | 2.3+ | GPL-2.0 |
| MMseqs2 | release_15 | GPL-3.0 |
| Snakemake | ≥ 8.0 | MIT |
| ChimeraX (visualisation) | ≥ 1.7 | UCSF non-commercial |
| Bio.PDB / Biopython | bundled | BSD |

## Reproducibility commands

The two canonical commands re-execute every quantitative claim in the manuscript.

```bash
# Tier 1 replay — no GPU required, < 30 seconds
python 06_Report/Mr_Repro/benchmark/run_replay.py

# RCSB co-crystal seed redocking — RTX 2060 (WSL2 Ubuntu 24.04), ~ 4 minutes
wsl -d Ubuntu -- bash -lc "cd /mnt/d/project/experiment/DAD && \
  python3 06_Report/Mr_Repro/benchmark/run_rcsb_redocking.py \
    --gnina 06_Report/Mr_Repro/tools/gnina/run_gnina_wsl.sh \
    --exhaustiveness 32 --num-modes 9 --seed 0"
```

Equivalent commands for Colab and Docker / Apptainer paths are listed in `06_Report/Mr_Pipeline/runtime_freeze.md` §2.

## Container image

`06_Report/Mr_Repro/container-hpc/Dockerfile` builds a GNINA-bundled image with pinned CUDA 12 / cuDNN 9 user-space libraries. Convertible to Apptainer SIF for HPC clusters.

## Continuous integration

`/.github/workflows/test.yml` runs the unit-test suite on every push (Tier 1 replay smoke test only; live GNINA tests are excluded from CI to avoid GPU dependency). pytest temp-directory permission notes for Windows hosts are in `runtime_freeze.md` §4.

## Conditional availability

The supporting primary paper (Sung J-Y et al., *iScience* 2026) is peer-reviewed and published, so the protocol's primary-paper requirement is already satisfied. The repository is currently in a private working tree pending submission of this protocol manuscript; the corresponding author has completed the GitHub OAuth step (2026-05-08) and the public repository (`dad-protocol`, proposed name) will be created and supplied to the editors at submission. At acceptance of this protocol manuscript the Zenodo archive will be minted.
