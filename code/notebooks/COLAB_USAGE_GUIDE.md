# DAD Protocol — Colab Usage Guide

**버전 / Version:** Phase E, 2026-05-08  
**대상 / Audience:** wet-lab biologists and computational biologists using Google Colab for the first time.

---

## §1 처음 5분 / First 5 Minutes

### 단계별 안내 / Step-by-step

| # | 한국어 | English |
|---|--------|---------|
| 1 | 배지 클릭 → 구글 계정으로 Colab 열기 | Click the "Open in Colab" badge → sign in |
| 2 | 메뉴 → **런타임 → 런타임 유형 변경 → GPU → T4** | Menu → **Runtime → Change runtime type → GPU → T4** |
| 3 | 메뉴 → **런타임 → 모두 실행** | Menu → **Runtime → Run all** |
| 4 | §0 Setup 셀이 condacolab을 설치하고 **커널을 자동 재시작**합니다 | §0 Setup installs condacolab; kernel restarts automatically |
| 5 | 재시작 완료 후 다시 **런타임 → 모두 실행** | After restart, click **Runtime → Run all** again |
| 6 | §1 Input 형식에서 기본값(Tier 1 replay) 확인 후 실행 완료 대기 | Check §1 Input defaults (Tier 1 replay) and wait for completion |

### GPU 활성화 확인 / Verify GPU

```python
import subprocess
result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
print(result.stdout[:200])
```

T4 또는 L4 GPU가 표시되면 정상입니다. GPU가 없으면 `exec_mode = "tier1_replay"` 로 유지하세요 (GPU 불필요).

---

## §2 Tier 1 Default vs Custom Input

### Tier 1 Default (권장 시작점 / Recommended starting point)

`USE_TIER1_DEFAULT = True` (기본값)로 설정하면:
- 단백질: NA23_RS01195 (MCP), NA23_RS08105 (CRP), NA23_RS00870 (RbsB) — 3개
- 리간드: Ala-Ile, Gly-Val — 2개 디펩타이드
- 케이스: 3 × 2 = 6 (Revision.txt ground truth 재현)
- GPU 불필요 (`exec_mode = "tier1_replay"`)

### Custom Input (본인 데이터)

`USE_TIER1_DEFAULT = False` 로 변경하면 §1 Input 폼의 custom 필드가 활성화됩니다.

**단백질 FASTA 형식:**
```
>ProteinA
MKTLLLSVALAGFASAHAA...
```
여러 단백질은 세미콜론(`;`)으로 구분:
```
>ProtA\nSEQ1;>ProtB\nSEQ2
```

**리간드 SMILES 형식 (name:SMILES;name:SMILES):**
```
Ala-Ile:CC(N)C(=O)NC(CC(C)C)C(O)=O;Gly-Val:NCC(=O)NC(C(C)C)C(O)=O
```

| 항목 | Tier 1 default | Custom |
|------|---------------|--------|
| 단백질 | 3 (AW-1 화학주성) | 1–20 권장 |
| 리간드 | 2 (디펩타이드) | 1–50 권장 |
| GPU 필요 | 아니오 | live 모드만 |
| 실행 시간 | ~2 분 | ~15 분/단백질 (ColabFold) + ~15 초/케이스 (GNINA) |

---

## §3 결과 해석 / Result Interpretation

### 수용 기준 / Acceptance threshold

| 지표 | 기준 | 해석 |
|------|------|------|
| top-pose RMSD | < 2.0 Å | PASS — 도킹 성공 |
| top-pose RMSD | ≥ 2.0 Å | FAIL — ranking_error 또는 search_failure 확인 |
| CNN pose score | 높을수록 좋음 | 전역 분류기 (AUROC 0.958) |
| vina score | 낮을수록 좋음 (kcal/mol) | 에너지 함수 |

### best-of-9 주의 / best-of-9 diagnostic note

> **best-of-9 RMSD는 탐색 능력 진단 지표입니다 — 운영 성공 지표가 아닙니다.**  
> best-of-9 is a **search-recovery diagnostic**, not an operational metric.  
> CNN score가 9개 포즈 중 최고 포즈를 1위로 선택하지 못하는 경우가 있습니다.  
> 이것이 12/16 (top-pose) vs 15/16 (best-of-9) gap의 원인입니다.

### 허용 오차 / Tolerance thresholds

GPU 아키텍처에 따라 GNINA 점수는 다음 범위 내에서 달라질 수 있습니다:

| 점수 | 허용 오차 |
|------|-----------|
| Vina score | ±0.5 kcal/mol |
| CNN pose score | ±0.05 |
| CNN affinity | ±0.5 pKi |

---

## §4 결과 다운로드 및 인용 / Download and Citation

### 결과 파일 / Output files

§12 Reproducibility 셀이 실행되면 자동 다운로드가 시작됩니다.

| 파일 | 내용 |
|------|------|
| `docking_results.tsv` | 모든 케이스의 점수 + RMSD |
| `table5_family_coverage.tsv` | Manuscript Table 5 |
| `table6_failure_cases.tsv` | Manuscript Table 6 |
| `fig1_rmsd_boxplot.{png,pdf}` | Figure 1 (300 dpi) |
| `fig2_cnn_rmsd_scatter.{png,pdf}` | Figure 2 (300 dpi) |
| `fig3_roc_curve.{png,pdf}` | Figure 3 (300 dpi) |
| `reproducibility_footprint.json` | Methods section 버전 정보 |
| `figures.zip` | 모든 그림 압축본 |

### 수동 다운로드 (자동 다운로드 실패 시)

```python
from google.colab import files
files.download('/content/DAD_run/results/docking_results.tsv')
```

### Methods 인용 문구 / Methods citation wording

```
The protein-metabolite docking analysis was performed using the DAD (Dynamic Affinity Dock)
protocol (Song and Lee, 2026), implemented as an interactive notebook hosted on GitHub and
executed on a cloud-based GPU environment via Google Colab, ensuring high-performance
computing without local hardware requirements. Protein structures were predicted using
ColabFold [Mirdita et al. 2022] and docking was performed with GNINA v1.3.2 [McNutt et al.
2021] with exhaustiveness=32, num_modes=9, seed=0.
```

---

## §5 문제 해결 / Troubleshooting

### libcudnn 오류 / libcudnn error

**증상:** `libcudnn.so.9: cannot open shared object file`  
**원인:** LD_LIBRARY_PATH 미설정  
**해결:** §0 Setup 셀이 자동으로 pip wheel cudnn 경로를 설정합니다. 셀을 다시 실행하세요.

```python
import os, site
site_pkgs = site.getsitepackages()[0]
os.environ['LD_LIBRARY_PATH'] = f'{site_pkgs}/nvidia/cudnn/lib'
```

### 세션 타임아웃 / Session timeout

**증상:** 연결 끊김 (90분 비활성 후)  
**해결:** 결과는 `/content/{job_name}/results/` 에 저장됩니다.  
Drive 마운트 후 `DRIVE_DAD_ROOT` 로 중간 결과를 저장하는 것을 권장합니다.

```python
# Drive 자동 저장 예시
import shutil
shutil.copytree('/content/DAD_run', '/content/drive/MyDrive/DAD/DAD_run', dirs_exist_ok=True)
```

### GPU 할당량 초과 / GPU quota exceeded

**증상:** "You've used all available compute units" 또는 T4 배정 불가  
**해결:** `exec_mode = "tier1_replay"` 로 설정하면 GPU 없이 실행 가능합니다. 또는 다음 날 재시도하세요.

### condacolab 재설치 오류 / condacolab reinstall error

**증상:** "cannot install condacolab in an already-initialized conda environment"  
**해결:** 메뉴 → **런타임 → 런타임 연결 해제 및 삭제** → 다시 연결 후 실행.

### GNINA 다운로드 실패 / GNINA download failure

**증상:** wget 오류 또는 binary 파일 크기 < 100 MB  
**해결:** GitHub releases에서 직접 다운로드하여 `/usr/local/bin/gnina` 에 업로드하세요.
URL: `https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2`

### Drive 마운트 후 경로를 찾지 못하는 경우

**증상:** `AW1_ref not found at /content/drive/MyDrive/DAD/AW1_ref`  
**확인:**
```python
from pathlib import Path
list(Path('/content/drive/MyDrive/DAD').iterdir())
```
`DRIVE_DAD_ROOT` 를 실제 경로에 맞게 수정하세요.

---

## §6 Methods 섹션 인용 / Methods Section Citation

### Computational Environment

```
The analysis was performed on a cloud-based GPU environment via Google Colab (Google LLC,
Mountain View, CA), ensuring high-performance computing without local hardware requirements.
```

### Software Versions

`reproducibility_footprint.json` 의 내용을 참조하여 Methods에 기재하세요:

```
Protein structure prediction: ColabFold v1.5.5 (Mirdita et al., Nat Methods 2022).
Molecular docking: GNINA v1.3.2 (McNutt et al., J Cheminform 2021), exhaustiveness=32,
num_modes=9, seed=0. Ligand 3D conformers: RDKit ETKDGv3 + MMFF optimization.
Pocket detection: P2Rank v2.4 (Krivak & Hoksza, J Cheminform 2018).
Python 3.10; pandas, numpy, scipy, scikit-learn (see footprint JSON for exact versions).
```

### Reproducibility Statement

```
All docking scores are reproducible within tolerance (vina ±0.5 kcal/mol, CNN pose ±0.05,
CNN affinity ±0.5 pKi) across GPU architectures. Bitwise identity is not assumed due to
GNINA's GPU floating-point non-determinism. The Tier 1 replay benchmark reproduces the
six values from the primary paper's reviewer-response docking table.
```
