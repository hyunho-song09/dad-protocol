# Two-Phase Architecture Plan — Structure Prep ↔ Selective Ligand Scoring

작성일: 2026-05-08 KST (Mode B post-iScience publication)
배경: 사용자 directive 2026-05-08 — AF2 가 Colab 에서 felt 너무 느려 protocol 의 사용성이 위협받음. 동시에 "FASTA 다수 + SMILES 다수" 패턴에서 매번 모든 조합을 dock 하는 것은 사용자 의도와 맞지 않음 (사용자는 protein subset × ligand subset 만 골라서 score 계산하기를 원함).

---

## 1. 결정 — Two-Phase

### Phase A — Structure Preparation (sequence-only, ligand 무관)
- 입력: multi-FASTA 파일 (header = unique protein name)
- 처리: 한 번의 `colabfold_batch` 호출 + 한 번의 P2Rank 일괄 실행
- 캐시: protein_name + sequence-SHA256 → reuse (재실행 시 재계산 X)
- 출력: `structure_registry.tsv` (protein_name / seq_sha / pdb_path / pocket_csv / pLDDT_mean / af2_runtime_s)

### Phase B — Selective Ligand Scoring (ligand 선택 + protein 부분집합 선택)
- 입력: multi-SMILES (with optional names) + ipywidgets 으로 2 multi-select (proteins, ligands)
- 처리: 선택된 protein × 선택된 ligand combination 만 GNINA 호출
- 캐시: case_id 별 결과 → 재실행 시 이미 계산된 pair skip
- 출력: `docking_master.csv` (선택 셀만), heatmap, ChimeraX cxc

**예시:**
- FASTA 입력: A, B, C, D, E (Phase A 한 번 실행 → 5 PDB + 5 pocket)
- SMILES 입력: F, G
- 사용자 선택 1: proteins={A,B,D}, ligands={F} → 3 dock 실행 (A-F, B-F, D-F)
- 사용자 선택 2: proteins={C}, ligands={F,G} → 2 dock (C-F, C-G), Phase A 재계산 X
- 사용자 선택 3: all × all → 10 dock (A-F, A-G, B-F, B-G, C-F, C-G, D-F, D-G, E-F, E-G)

---

## 2. AF2 Speed Strategy

`AF2_colabfold_batch_setting_plan.md` (사용자 제공) 의 권장값 기반 + 추가 기법:

| 기법 | 설정 | 효과 |
|------|------|------|
| Single rank model | `--num-models 1` | **~5× 빠름** vs 5-model default |
| Reduced recycles | `--num-recycle 3` | 기본 3 유지 (1 로 더 빠르나 정확도↓) |
| Length-sorted batch | `--sort-queries-by length` | JAX/XLA compilation cache 재사용 (큰 효과) |
| MMseqs2 API MSA | `--msa-mode mmseqs2_uniref_env` | 로컬 MSA 빌드 회피 |
| Skip AMBER relax | `--num-relax 0` (option) | 추가 ~2× 빠름 (사용자 선택) |
| Templates off | `--templates 0` (option) | template fetch I/O 절약 |
| Batch JAX cache | (자동, batch mode 의 효과) | 1 query 당 압축률 유지 |
| Single command | 5 query 1 run | colabfold_batch 의 핵심 |

**Conservative default (DAD Phase A):** `--num-models 1 --num-recycle 3 --msa-mode mmseqs2_uniref_env --sort-queries-by length --zip`

**Fast preview option (사용자 toggle):** 위 + `--num-relax 0 --templates 0`

Multi-mer 입력은 `:` chain break 그대로 지원 (`--model-type alphafold2_multimer_v3` 자동 감지).

---

## 3. Cache Architecture

### Phase A cache directory
```
WORK_DIR/
├── phase_a/
│   ├── input_sequences.fasta              # validated multi-FASTA
│   ├── colabfold_out/                     # raw colabfold_batch output
│   │   ├── <name1>_unrelaxed_rank_001*.pdb
│   │   └── ...
│   ├── structure_cache/                   # canonical PDBs keyed by SHA
│   │   ├── <seq_sha[:12]>_<name>.pdb
│   │   └── ...
│   ├── pocket_cache/                      # P2Rank outputs keyed by SHA
│   │   ├── <seq_sha[:12]>_<name>_predictions.csv
│   │   └── ...
│   └── structure_registry.tsv             # name | seq_sha | pdb | pocket | pLDDT | runtime
```

### Phase B cache
```
WORK_DIR/
└── phase_b/
    ├── ligands/                           # 3D SDFs keyed by SMILES SHA
    ├── runs/                              # per-pair docking runs keyed by case_id
    │   ├── <protein>_<ligand>/
    │   │   ├── docked.sdf
    │   │   ├── scores.json
    │   │   └── interactions.tsv
    └── docking_master.csv                 # selected pairs only, append-only
```

### Hash key
- Protein cache key: `sha256(sequence_string)[:12]`
- Ligand cache key: `sha256(canonical_smiles)[:12]`
- Pair cache key: `<protein_sha>_<ligand_sha>` (case_id)

Re-run safety: 새 SMILES 추가 시 Phase A 0% 재계산. 새 protein 추가 시 새 protein 만 계산.

---

## 4. UI Pattern (ipywidgets)

Phase B 의 selection cell:

```python
import ipywidgets as widgets
from IPython.display import display

protein_select = widgets.SelectMultiple(
    options=list(structure_registry.index),
    description='Proteins:',
    rows=8,
)
ligand_select = widgets.SelectMultiple(
    options=list(ligand_registry.index),
    description='Ligands:',
    rows=8,
)
all_proteins_btn = widgets.Button(description='Select all proteins')
all_ligands_btn = widgets.Button(description='Select all ligands')
run_btn = widgets.Button(description='Run docking on selection', button_style='primary')

display(protein_select, ligand_select, all_proteins_btn, all_ligands_btn, run_btn)

# button callbacks: select all, run combinatorial dispatcher
```

---

## 5. AF3 Path Coexistence

AF3 mode (사용자 directive 2026-05-08 retained as option):
- `STRUCTURE_MODE = "af3_results"` 가 설정되면 `AF3_RESULTS_PATH` 의 .cif → `structure_registry.tsv` 로 직접 ingestion (Phase A 의 ColabFold 호출 우회)
- 그 외 모든 STRUCTURE_MODE 옵션은 Phase A 의 ColabFold/ESMFold 경로로 통합

User-facing 결과: AF3 가 Colab 에서 직접 실행은 부적합 (DeepMind 의 restricted weights), 사용자 외부 (alphafoldserver.com 또는 로컬 install) 결과를 받아오는 path 만 유지.

---

## 6. STRUCTURE_MODE 갱신 (5 → 5 옵션, default 변경)

| 모드 | 용도 | 변경점 |
|------|------|--------|
| `"colabfold_af2"` | **(NEW DEFAULT)** Phase A 의 colabfold_batch 호출 | 사용자 directive 정합 |
| `"af3_results"` | 사용자가 alphafoldserver.com 결과 보유 시 | (유지) |
| `"esmfold_api"` | 짧은 protein API fallback | (유지) |
| `"user_pdb"` | 사용자 직접 PDB 업로드 | (유지) |
| `"auto"` | cache → ESMFold → ColabFold AF2 chain | (유지, 코덱스 fallback) |

Default 변경 이유: AF3 직접 실행 불가 확인 (사용자 2026-05-08 보고). Phase A 가 ColabFold 중심으로 동작하는 게 자연스러움.

---

## 7. 코덱스 fix 보존 (7 commit, 1317253 ~ 7caa1df)

이번 재설계는 코덱스의 다음 fix 를 모두 그대로 유지:
- condacolab restart `RESTART_REQUESTED` flag (ae23ba9)
- `colabfold[alphafold]` extras 자동 install (06096ff)
- Python 3.12 ColabFold TensorFlow `.so` cleanup 패턴 (5c78361)
- ColabFold AF2 fallback 의 install 명령 (`colabfold[alphafold-minus-jax] @ git+...`) (5c78361)
- AF2 relax dependencies (OpenMM/PDBFixer) install (7caa1df)
- ESMFold short-protein API fallback (5c78361)
- user-input focused workflow 패턴 (f66bca5)

새로 추가:
- `colabfold_batch` CLI mode (multi-FASTA 입력) — single-query mode 와 공존
- length-sorted batch 옵션
- `--num-relax 0 --templates 0` 의 fast-preview toggle
- structure_registry.tsv + sequence-SHA cache
- ipywidgets selection UI
- Phase A / Phase B 명확한 셀 분할

---

## 8. Notebook Cell Layout (재설계)

### Phase A cells (FASTA → PDB + Pocket)
- §0 Setup environment (코덱스 작업 보존)
- §1 Phase A input — multi-FASTA Colab Form
- §2 Phase A run — `colabfold_batch` (또는 AF3 ingestion / ESMFold / user_pdb)
- §3 Phase A pocket detection — P2Rank batch
- §4 Phase A summary — structure_registry.tsv + per-protein pLDDT bar

### Phase B cells (SMILES → selective dock)
- §5 Phase B input — multi-SMILES Colab Form (`name:SMILES;...`)
- §6 Phase B selection UI — ipywidgets multi-select (proteins × ligands)
- §7 Phase B run — combinatorial GNINA dispatcher
- §8 Phase B aggregate — docking_master.csv + selected-pair heatmap
- §9 Phase B visualize — py3Dmol + ChimeraX cxc
- §10 Reproducibility footprint — manifest.json (Phase A + Phase B SHA256)

---

## 9. Manuscript Reflection (Mr_Prof 작업)

§7 Procedure 를 두 개로 분리:
- §7A Phase A — Structure preparation (Stages 0–5)
- §7B Phase B — Selective ligand scoring (Stages 6–12, but with selection step)

§6 Equipment setup 의 Stage 4 row 갱신:
- AF2/ColabFold default + `--num-models 1` rationale + `colabfold_batch` multi-FASTA
- AF3 alternative (사용자 directive 변경 반영)
- ESMFold short-protein fallback

§9 Anticipated Results §10 의 worked example 도 두 phase 로 나눠 시연 (Tier 1 6 case 의 결과는 동일).

---

## 10. Codex Tracking 포인트

코덱스가 GitHub 에서 직접 commit 추가할 때 참고할 수 있도록 다음 항목을 체크리스트화:

| 영역 | 확인 |
|------|------|
| Phase A `colabfold_batch` CLI 호출 정상 | `--num-models 1 --num-recycle 3 --sort-queries-by length --msa-mode mmseqs2_uniref_env --zip` |
| sequence SHA cache hit/miss 로직 | `WORK_DIR/phase_a/structure_cache/` 의 파일 이름 패턴 |
| P2Rank batch 호출 | structure_registry 의 모든 protein 에 대해 한 번씩 |
| ipywidgets 의존성 | Colab 기본 설치되어 있어야 함 (확인) |
| Phase B 의 case_id 일관성 | `<protein_sha[:12]>_<ligand_sha[:12]>` |
| AF3 mode coexistence | `STRUCTURE_MODE='af3_results'` 일 때 §2 cell skip |
| 코덱스 7 fix 보존 | 위 §7 의 7 항목 |
| docking_master.csv append-only | concurrent run 시 conflict 방지 |
| sub-tasks 의 progress bar | tqdm 사용 (사용자 directive: progress bar 추가됨) |

---

## 11. 다음 단계 — 작업 분배

| 에이전트 | 작업 |
|---------|------|
| Mr_Pipeline | 본 계획서 따라 notebook (top-level + code/notebooks/) 재구성. 2-phase 셀 분할 + colabfold_batch + selection UI + cache 로직 |
| Mr_Prof | manuscript §6 / §7 / §9 의 two-phase reflection (소규모) |
| Mr_Repro | (옵션) WORK_DIR layout 의 cache 일관성 verification 스크립트 |
| Mr_Dock | (옵션) GNINA dispatcher 의 dispatch 로직 검증 (이미 docking.py 구현 있음) |
| team-lead | commit + push + 사용자 검증 follow-up |

본 계획서는 codex 가 GitHub `development_logs/two_phase_architecture_plan.md` 에서 직접 읽고 추적 가능. Mr_Pipeline 의 notebook 변경 사항은 모두 본 plan 의 섹션 번호로 cross-reference 됨.
