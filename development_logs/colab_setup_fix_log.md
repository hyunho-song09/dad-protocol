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
