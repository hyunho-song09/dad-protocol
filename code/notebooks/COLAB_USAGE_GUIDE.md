# DAD Colab Usage Guide

## Quick Start

1. Open `DAD_protocol.ipynb` in Colab.
2. Use GPU only for `live` mode.
3. Click `Runtime > Run all`.
4. If Setup installs condacolab, the kernel restarts.
5. After restart, click `Runtime > Run all` again.

## Recommended First Test

Use the defaults:

- `USE_TIER1_DEFAULT = True`
- `exec_mode = "tier1_replay"`
- `INSTALL_GNINA = True` or `False`

Tier 1 replay does not need GPU or live GNINA execution.

## Setup Cell Behavior

The Setup cell does four things:

- installs condacolab on first Colab run;
- stops after the required kernel restart;
- installs RDKit/Biopython after restart;
- downloads and checks GNINA only when `INSTALL_GNINA=True`.

Expected success messages:

- `Package check: PASS`
- `GNINA check: PASS` when GNINA is installed

## Input Modes

| Mode | Use case | GPU |
|---|---|---|
| `tier1_replay` | Fast smoke test and deterministic replay | No |
| `live` | Real GNINA docking and structure workflow | Yes |

## Custom Input

Protein FASTA:

```text
>ProteinA
MKTLLLSVALAGFASAHAA...
```

Multiple proteins:

```text
>ProtA\nSEQ1;>ProtB\nSEQ2
```

Ligand SMILES:

```text
Ala-Ile:CC(N)C(=O)NC(CC(C)C)C(O)=O;Gly-Val:NCC(=O)NC(C(C)C)C(O)=O
```

## Main Outputs

| File | Content |
|---|---|
| `docking_results.tsv` | scores and RMSD by case |
| `table5_family_coverage.tsv` | family summary |
| `table6_failure_cases.tsv` | failure-case summary |
| `figures.zip` | generated figures |
| `reproducibility_footprint.json` | runtime and package metadata |

## Troubleshooting

### Setup stops after condacolab

Expected. The kernel was restarted. Run all cells again.

### mamba install fails

The notebook now stops before mamba is called in the old kernel. After restart, rerun all cells. If mamba still fails, the cell falls back to pip for RDKit/Biopython.

### `libcudnn.so.9` error

Rerun Setup with `INSTALL_GNINA=True`. The cell installs CUDA/cuDNN wheels and sets `LD_LIBRARY_PATH`.

### GPU quota exceeded

Use:

```python
exec_mode = "tier1_replay"
```

### GNINA download fails

Rerun Setup. If it still fails, manually place the binary at:

```text
/usr/local/bin/gnina
```

Release URL:

```text
https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2
```

## Interpretation

- Top-pose RMSD < 2.0 A: PASS.
- Best-of-9 RMSD is diagnostic only.
- CNN pose score is useful across cases, but can rank a non-best pose first within one case.

## Tolerances

| Score | Tolerance |
|---|---|
| Vina | +/-0.5 kcal/mol |
| CNN pose | +/-0.05 |
| CNN affinity | +/-0.5 pKi |
