# Figure 3 References

## Primary figure code references used

### 1. PoseBusters — Chemical Science 2024 (primary boxplot reference)
- **Paper**: Buttenschoen A, Morris GM, Deane CM. PoseBusters: AI-based docking methods fail to generate physically valid poses or generalise to novel sequences. *Chem Sci* 2024;15:3130–3139. doi:10.1039/D3SC04185A
- **GitHub**: https://github.com/maabuu/posebusters
- **Applied**: grouped boxplot layout with family grouping, 2.0 Å threshold line, FAIL case annotation with red dots

### 2. DiffDock — ICLR 2023 (RMSD percentile reference)
- **Paper**: Corso G, Stark H, Jing B, Barzilay R, Jaakkola T. DiffDock: Diffusion steps, twists, and turns for molecular docking. *ICLR* 2023.
- **GitHub**: https://github.com/gcorso/DiffDock
- **Applied**: RMSD threshold annotation style, top-pose vs best-of-9 dual reporting convention

## Data source
- `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_results.tsv` — 16 rows
- `06_Report/Mr_Repro/external_validation/rcsb_seed/redocking_summary.json` — family aggregate

## Style decisions
- Top-pose: blue (#0072B2), best-of-9: sky blue (#56B4E9)
- Failure annotation: vermillion (#D55E00) triangles + arrows
- Threshold: red dashed (#CC0000)
- Jittered dots overlaid on boxplots for individual case visibility
- n= annotation per family group
