# Colab Setup Fix Log

## 2026-05-08: condacolab restart shown as SystemExit traceback

### User-facing symptom

In Google Colab, `0. Setup Environment` installed `condacolab`, restarted the kernel, and then displayed:

```text
SystemExit: Kernel restart requested. After restart, click Runtime > Run all.
```

This looked like a notebook error even though the restart itself was expected.

### Root cause

The setup cell intentionally raised `SystemExit` after `condacolab.install()` to prevent the old pre-Conda kernel from continuing into `mamba` package installation. Colab displays `SystemExit` as a traceback-like exception.

### Fix

The setup cell now uses a `RESTART_REQUESTED` flag:

- install `condacolab`;
- call `condacolab.install()`;
- mark setup progress as `restart requested`;
- print a normal instruction message;
- skip the rest of the setup cell without raising `SystemExit`.

After Colab reconnects, the user should run all cells again. On the second run, `mamba` is available and the rest of setup proceeds.

### Files updated

- `DAD_protocol.ipynb`
- `code/notebooks/DAD_protocol.ipynb`

### Related previous fixes

- The notebook was changed from paper-replay-first behavior to user-input-first behavior.
- Raw protein sequences without FASTA headers are accepted as `Protein_1`.
- Raw SMILES without `name:` prefixes are accepted as `Ligand_1`.
- Step-level progress bars were added throughout the notebook.

## 2026-05-08: ColabFold installed without AlphaFold extras

### User-facing symptom

`3. Structure Input or Prediction` failed in Colab with:

```text
ModuleNotFoundError: No module named 'alphafold'
RuntimeError: alphafold is not installed. Please run `pip install colabfold[alphafold]`
```

### Root cause

The notebook only checked whether the `colabfold_batch` executable existed. In the failing Colab runtime, `colabfold_batch` existed but the optional `alphafold` dependency was missing.

### Fix

The setup and structure cells now treat ColabFold as ready only when both conditions are true:

- `colabfold_batch` is available;
- the Python module `alphafold` is importable.

If either condition fails, the notebook installs:

```bash
pip install -q "colabfold[alphafold]"
```

The structure cell writes installation diagnostics to:

```text
WORK_DIR/structures/colabfold_install.log
```

and ColabFold run diagnostics to:

```text
WORK_DIR/structures/colabfold.log
```


## 2026-05-08: ColabFold TensorFlow binary import failure on Python 3.12

### User-facing symptom

`3. Structure Input or Prediction` failed after installing AlphaFold extras:

```text
ImportError: ... _pywrap_tensorflow_lite_metrics_wrapper.so: undefined symbol
RuntimeError: ColabFold failed. See structures/colabfold.log
```

### Root cause

The failing runtime reached AlphaFold/TensorFlow import, but the TensorFlow binary loaded by Colab Python 3.12 was incompatible. This is not a GPU quota error and is not fixed by repeatedly reinstalling `colabfold[alphafold]`.

### Fix

- Default structure mode is now `esmfold_api`, which posts the input sequence to the ESMFold PDB API and avoids AlphaFold/TensorFlow installation.
- `colabfold` remains an optional mode, but the notebook blocks it on Python 3.12+ with an explicit error message.
- Predicted or uploaded PDB files are cached by protein sequence hash.
- Users can keep the same `job_name`, change `custom_ligand_smiles`, and rerun from `1. Input Data` onward to score new ligands against the cached PDB.

### Files updated

- `DAD_protocol.ipynb`
- `code/notebooks/DAD_protocol.ipynb`
- `README.md`
- `code/notebooks/COLAB_USAGE_GUIDE.md`
