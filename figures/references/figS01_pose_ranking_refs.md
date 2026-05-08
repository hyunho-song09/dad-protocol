# Supplementary Figure S1 References

## Primary figure code references used

### 1. PoseBusters — Chemical Science 2024 (supplementary bar plot reference)
- **Paper**: Buttenschoen A, Morris GM, Deane CM. PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. *Chem Sci* 2024;15:3130–3139. doi:10.1039/D3SC04185A
- **GitHub**: https://github.com/maabuu/posebusters
- **Applied**: case-level bar plot with categorical x-axis, binary color coding (PASS/FAIL), hatch pattern for secondary annotation (ranking gap)

## Codex Round 6 language constraint (C-R6-3)
- best-of-9 RMSD labeled "diagnostic" only, not "operational"
- wording: "within-case ranking can place a bad pose first even when good pose exists"

## Data source
- `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_results.tsv`
- Columns: case_id, best_pose_index, redock_pass, top_pose_rmsd_a, best_pose_rmsd_a

## Key finding documented
- 4/16 cases (25%) have best_pose_index > 1:
  - TAR_ASP_2LIG: index 2 (top RMSD=0.797Å, best RMSD=0.354Å)
  - LBP_PHE_1USI: index 2 (top RMSD=0.428Å, best RMSD=0.404Å)
  - LBP_LEU_1USK: index 2 (top RMSD=0.799Å, best RMSD=0.463Å)
  - RBSB_RIB_3KSM: index 4 (top RMSD=12.778Å, best RMSD=1.572Å) — ranking failure

## Style decisions
- PASS bars: green (#009E73), FAIL bars: vermillion (#D55E00)
- Ranking gap hatching: /// pattern on bars with index > 1
- Annotation: per-bar RMSD values, ranking gap arrows for index > 1 cases
