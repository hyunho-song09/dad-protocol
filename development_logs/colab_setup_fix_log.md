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


## 2026-05-08: ESMFold 413 and ColabFold Python 3.12 fallback

### User-facing symptom

Longer custom proteins failed with:

```text
HTTPError: HTTP Error 413: Request Entity Too Large
```

Switching to `colabfold` also failed because the previous notebook blocked ColabFold on Python 3.12.

### Root cause

The ESMFold public API can reject long requests. The previous ColabFold block was too conservative: the project reference `AW1_ref/AlphaFold2.ipynb` already includes the ColabFold v1.6.1 install path and a TensorFlow shared-object cleanup that avoids the observed TensorFlow import crash.

### Fix

- `STRUCTURE_MODE` now defaults to `auto`.
- Auto mode tries cache/uploaded PDB first, then ESMFold API for shorter proteins, then ColabFold AF2 fallback.
- ColabFold install now uses the AlphaFold2.ipynb backend: `colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold`.
- The notebook applies the TensorFlow `.so` cleanup from the reference AlphaFold2 notebook before ColabFold import.
- AF2 quality controls were added: `AF2_PRESET`, `AF2_MODEL_TYPE`, `AF2_NUM_RECYCLES`, `AF2_NUM_MODELS`, `AF2_NUM_SEEDS`, and `AF2_RELAX_TOP_N`.
- If `AF2_RELAX_TOP_N > 0`, the structure step installs OpenMM/PDBFixer via mamba before AMBER relaxation.


## 2026-05-08: AF3-default + AF2-option structure prediction

### User-facing change

`STRUCTURE_MODE` now defaults to `"af3_results"` instead of `"auto"`. The notebook now has two structure cells:

- **§3a (new default):** Loads AlphaFold 3 CIF results from a user-specified folder (`AF3_RESULTS_PATH`) or a Colab ZIP upload. Parses per-model confidence JSON (ipTM, pTM, mean pLDDT, ranking_score), converts CIF → PDB via Bio.PDB / gemmi fallback, and renders a per-residue pLDDT bar chart with 70 / 90 threshold lines.
- **§3b (optional):** Retains the full AF2 / ESMFold / user_pdb path from previous codex commits. Skipped automatically when all structures are already loaded by §3a.

### STRUCTURE_MODE options (5)

| Value | Meaning |
|---|---|
| `af3_results` | Load CIF results from AF3 Server or local AF3 install (new default) |
| `auto` | Try cache/upload → ESMFold API → ColabFold AF2 |
| `user_pdb` | Upload or paste a PDB file directly |
| `esmfold_api` | ESMFold public API (≤ 400 aa) |
| `colabfold_af2` | ColabFold AF2 in-notebook prediction |

### New Colab Form params (§3a)

- `AF3_RESULTS_PATH` — path to the AF3 job folder (string, empty = Colab upload prompt)
- `AF3_USE_TOP_RANKED` — `True` uses model_0 only; `False` loads all 5 models and picks highest ranking_score

### Codex fixes preserved

All 7 prior codex fixes were preserved intact in §3b: RESTART_REQUESTED condacolab pattern, TF `.so` cleanup (`patch_tensorflow_for_colabfold()`), `ensure_colabfold_ready()` with `colabfold[alphafold-minus-jax]`, `apply_af2_preset()`, `ensure_amber_ready()`, `AF2_RELAX_TOP_N`, and the ESMFold 413 fallback chain.

### Files updated

- `Publication/DAD_protocol.ipynb`
- `Publication/code/notebooks/DAD_protocol.ipynb`
- `Publication/code/README.md`
