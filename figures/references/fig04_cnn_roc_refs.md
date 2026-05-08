# Figure 4 References

## Primary figure code references used

### 1. GNINA — J. Cheminformatics 2021 (score calibration scatter reference)
- **Paper**: McNutt AT, Francoeur P, Aggarwal R, Masuda T, Meli R, Ragoza M, Sunseri J, Koes DR. GNINA 1.0: molecular docking with deep learning. *J Cheminform* 2021;13:43. doi:10.1186/s13321-021-00522-2
- **GitHub**: https://github.com/gnina/gnina
- **Applied**: CNN score vs RMSD scatter plot layout, two-panel figure (a) scatter (b) ROC

## Quantitative values (Codex Round 6 — authoritative)
- Pearson r (top-pose RMSD vs CNN pose score): −0.889
- Pearson r (best-of-9 RMSD vs CNN pose score): −0.388
- AUROC: 0.958
- AUPRC: 0.987
- EF@25%: 1.333

## Data source
- `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_results.tsv`
- Columns used: cnn_pose_score, top_pose_rmsd_a, best_pose_rmsd_a, redock_pass

## Style decisions
- Scatter: PASS=green, FAIL=vermillion triangle
- Regression line: dashed gray
- ROC: blue curve, gray dashed reference diagonal
- AUROC text box inset (lower-right)
- 2-panel layout: 1.2:1.0 width ratio
