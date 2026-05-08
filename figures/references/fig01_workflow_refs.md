# Figure 1 References

## Primary figure code references used

### 1. ColabFold — Nature Methods 2022
- **Paper**: Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M. ColabFold: making protein folding accessible to all. *Nat Methods* 2022;19:679–682. doi:10.1038/s41592-022-01488-1
- **GitHub**: https://github.com/sokrypton/ColabFold
- **Figure code inspected**: `ColabFold/utils/plot_scores.ipynb` — pLDDT/PAE plot style, matplotlib rcParams, color conventions
- **Applied**: horizontal stage-box layout with phase color bands modelled on ColabFold Fig 1 workflow style

### 2. AlphaFill — Nature Methods 2023
- **Paper**: Hekkelman ML et al. AlphaFill: enriching AlphaFold models with ligand cofactors and ions. *Nat Methods* 2023;20:205–213. doi:10.1038/s41592-022-01685-y
- **GitHub**: https://github.com/PDB-REDO/alphafill
- **Applied**: stage-box + decision diamond layout for pipeline schematic

### 3. Nature Protocols figure standards
- **Reference**: Nature Portfolio author guidelines, 2025. https://www.nature.com/nprot/for-authors/preparing-your-submission
- **Applied**: Arial font, 183 mm double-column width, dpi=300, colorblind-safe palette (Wong 2011)

## Style decisions
- Phase palette: Wong 2011 colorblind-safe 8-color (Blue=#0072B2, Vermillion=#D55E00, Green=#009E73, Pink=#CC79A7, SkyBlue=#56B4E9, Yellow=#F0E442)
- Font: Arial 7pt body, 5.5pt annotations
- Layout: horizontal left-to-right, 4 phase bands, decision diamond at Stage 3
- DeepTMHMM EXCLUDE branch shown below Stage 3 decision diamond
