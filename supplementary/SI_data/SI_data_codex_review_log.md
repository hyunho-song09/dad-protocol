# SI Data — Codex External Review Log (Audit Trail)

본 문서는 DAD 의 외부 codex 리뷰 6 라운드와 DAD 팀 응답의 timeline 요약. 감사 추적용 SI.

## Round 1 — Initial codex review (2026-05-07)
- 입력: codex-feedback/strict_feedback.md + priority_action_plan.md
- 10 critical findings (C1–C10) + 7 milestones (M0–M7) + 7 strengths (S1–S7)
- 핵심 지적: end-to-end pipeline 미작동 / core ↔ Snakemake 미연결 / TriageDecision schema mismatch / Tier 1 재현 결과 부재 / AW1_ref 미관리 / env 무거움 / tests not clean / manuscript 앞서감 / repo not publication-ready
- DAD 응답: 5 에이전트 병렬 (Mr_Bio C3 / Mr_Repro C4-C7-C10 + M0-M1-M6 / Mr_Prof C9 / Mr_Struct C6 / Mr_Dock C5)

## Round 2 — Mr_Pipeline integration verified (2026-05-07)
- DAD 응답: io.py TriageDecision.triage_status optional field / Snakefile NotImplementedError 0 hits / replay mode `execution.mode=tier1_replay` / aw1_ref manifest-driven
- team-lead 독립 grep 검증 통과

## Round 2 follow-up (codex_round2_followup_report.md, 2026-05-07)
- codex 직접 4 auto-fix: utf-8-sig BOM (run_replay + structure.load_aw1_asset), Snakefile Revision.txt config-key, config.yaml benchmark.revision_table_path
- DAD 응답: 4 auto-fix 채택 + Phase D scope-out (data/ingest_phase_d/ 8 modules + benchmark/scope_phase_d/run_all.skeleton) + codex_response_round2.md

## Round 3 — RCSB seed external validation (2026-05-08 KST)
- codex 직접: build_rcsb_seed_validation.py + 16 case bridge dataset + family-balanced selection + contact_quality 16/16 PASS + redocking_queue.tsv

## Round 4 — Redocking preparation (2026-05-08 KST)
- codex 직접: prepare_rcsb_redocking.py + ligand atom completeness fix (1MDP→3MBP) + GNINA 1.3.2 binary download + 16/16 prep PASS
- 이 시점 blocker: WSL libcudnn.so.9 missing

## Round 5 — Live GNINA redocking (2026-05-08 KST)
- codex 직접: WSL CUDA/cuDNN runtime (tools/gnina_runtime/) + run_gnina_wsl.sh wrapper + GNINA score parser fix (CNNscore tag vs value) + RMSD evaluator fix (Kabsch removal → in-place RDKit symmetry-corrected) + run_rcsb_redocking.py --case-id/--limit
- 결과: 16/16 GNINA execution PASS, 12/16 top-pose RMSD ≤ 2.0 A, 15/16 best-of-9 RMSD ≤ 2.0 A, median top 0.561 A, median best-of-9 0.385 A

## Round 6 — Nature Protocols readiness review (2026-05-08 KST)
- codex 직접 7 auto-fix: primary_paper_relationship.md "satisfied" → "conditional" / manuscript Tier 2 wording → "planned Phase D" / failure_analysis.md best-of-9 → diagnostic only / runtime_freeze.md Ubuntu 22.04 → 24.04.3 LTS / gnina_runtime/MANIFEST.md added / sugar_sbp_extension_candidates.md current-status table corrected (7 actual cases) / enrichment_metrics.py 5 functions actually implemented (smoke test pass)
- codex 정량 분석: CNN pose vs top RMSD Pearson r=-0.889, Spearman ρ=-0.656; CNN as PASS/FAIL classifier AUROC=0.958, AUPRC=0.987, EF@25%=1.333
- 핵심 framing: "within-case ranking can place a bad pose first even when good pose exists in generated 9 modes"
- 5 Blockers: (B1) primary paper DOI / (B2) one live end-to-end / (B3) clean reproducibility log / (B4) reference-free pose ranking / (B5) validation scale
- DAD 응답: codex_response_round6.md + 4 에이전트 병렬 (Mr_Pipeline DAD_protocol.ipynb / Mr_Dock select_consensus_pose 5-method / Mr_Repro sugar SBP 5 prep + WSL snakemake plan + enrichment plan / Mr_Prof manuscript R6 reflection + presubmission_enquiry_summary.md)

## Cumulative status (post-R6)
- Critical findings closed (verified): C1, C2, C3, C5, C6, M2, M3, M4
- Closed (unverified): C4, C7, C9, C10, M0, M1, M6, M7
- Deferred to user: C8 (env-lite pytest)
- Phase D in progress: B2 user execution / B3 WSL snakemake / B5 sugar SBP 5

## Files in codex-response/
- codex_round{2_followup,3_external_validation,4_redocking_preparation,5_redocking_runtime_resolution,6_claude_review_nature_protocols_readiness}_*.md (codex 자체 보고)
- codex_round{2,3,4,5}_claude_handoff.md (codex 가 Claude 에게 다음 작업 명세)
- codex_round2_test_log.md (codex 자체 테스트 로그)
- codex_response_round2.md, codex_response_round5.md, codex_response_round6.md (DAD 팀 응답)
- adoption_log.md (finding × 응답 매트릭스, append-only)
- artifacts_index.md (finding → 변경 파일 deep-link)
- verification_results.md (객관적 실행 증거 R1-R4)
- next_review_request.md (다음 codex 리뷰 ask)
- claude_access_index.md (Claude 가 우선 읽을 경로 + ignore 영역)
