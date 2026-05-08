# DAD Test Suite

**Location:** `code/tests/`  
**Total tests:** ~134-135 depending on optional dependency skips (33 triage + 56 docking/interaction + 45 structure/pocket plus collection-level helpers)

---

## Running Tests

### Prerequisites

```bash
conda create -n dad-lite -c conda-forge python=3.11 rdkit biopython numpy scipy pandas pytest
conda activate dad-lite
```

No GPU, GNINA binary, ColabFold, P2Rank, or network access required.
All external tool calls are mocked.

### Run all tests

```bash
cd Publication/code
pytest tests/ -v
```

### Run individual test files

```bash
# Stage 3 triage tests (33 tests)
pytest tests/test_triage.py -v

# Stage 6-10 docking/interaction/visualization tests (56 tests)
pytest tests/test_docking.py -v

# Stage 4-5 structure/pocket tests (45 tests)
pytest tests/test_structure.py -v
```

---

## Test Coverage

### `test_triage.py` — 33 tests

Tests for `dad.core.triage` (Stage 3 pre-docking biological triage):

| Group | Tests | What is validated |
|-------|-------|-------------------|
| R1 length filter | 3 | EXCLUDE < 50 aa, PASS >= 50 aa |
| R2 signal peptide | 4 | SP clipping, mature sequence coords |
| R3 TM topology | 8 | nTM=0 PASS, nTM=2 domain extraction, nTM≥7 EXCLUDE |
| R4 dock region | 3 | Minimum dock-region length enforcement |
| R5 functional class | 5 | Pfam-to-class mapping, priority assignment |
| TriageRecord export | 4 | TSV writer, FASTA writer, verdict mapping |
| Tier 1 proteins | 6 | MCP (nTM=2, SP clip), Crp (nTM=0), RbsB (SP+SBP) |

Mock topologies derived from Phobius/DeepTMHMM literature values.  
No network calls are made.

### `test_docking.py` — 56 tests

Tests for `dad.core.docking`, `dad.core.interaction`, `dad.core.visualize`:

| Group | Tests | What is validated |
|-------|-------|-------------------|
| Ligand prep | 8 | SMILES→SDF, RDKit fallback, max_dim computation |
| Box sizing | 6 | DAD box rule: max(22, max_dim+10) |
| GNINA SDF parsing | 6 | Multi-model SDF, REMARK score extraction |
| Replay mode | 6 | `replay_from_ground_truth()`, alias normalisation |
| Consensus pose | 6 | ensemble/cluster/posebusters/plip/consensus methods |
| Interaction profiling | 8 | Contact extraction, z-score ranking, master CSV |
| Visualization | 6 | HTML generation, ChimeraX .cxc script |
| Integration | 10 | End-to-end io schema round-trips |

### `test_structure.py` — 45 tests

Tests for `dad.core.structure` and `dad.core.pocket`:

| Group | Tests | What is validated |
|-------|-------|-------------------|
| pLDDT parsing | 6 | B-factor column → mean pLDDT, low_confidence_flag |
| AW1_ref reuse | 8 | `load_existing_pdb()`, `load_aw1_asset()` manifest resolution |
| AF3 support | 6 | CIF parsing, multi-model loading, pLDDT from JSON |
| Dock region trim | 5 | Bio.PDB select + plain-text fallback |
| P2Rank parsing | 8 | AW1_ref CSV format, space-padded fields |
| Pocket-to-box | 6 | Box size rule, negative ligand_max_dim guard |
| Integration | 6 | StructurePrediction → PocketResult → DockingBox chain |

---

## Expected Results

All tests pass on:
- Python 3.10+
- rdkit ≥ 2023.09
- numpy ≥ 1.24
- biopython ≥ 1.81

Approximate runtime: < 30 seconds on any modern CPU.

---

## Notes

- Tests use `unittest.mock` to patch GNINA, ColabFold, P2Rank subprocess calls.
- AW1_ref file paths in tests are resolved relative to the project root;
  tests requiring AW1_ref files are skipped with `pytest.mark.skipif` when
  running outside the project directory.
- `conftest.py` adds `Publication/code/` to `sys.path` automatically.
