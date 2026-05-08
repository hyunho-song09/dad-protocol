# Figure 5 References

## Primary figure code references used

### 1. ColabFold — Nature Methods 2022 (protocol figure style)
- **Paper**: Mirdita M et al. ColabFold: making protein folding accessible to all. *Nat Methods* 2022;19:679–682. doi:10.1038/s41592-022-01488-1
- **GitHub**: https://github.com/sokrypton/ColabFold
- **Applied**: rainbow N→C cartoon coloring convention, compact figure layout with score annotation box

### 2. AlphaFill — Nature Methods 2023 (interaction view style)
- **Paper**: Hekkelman ML et al. AlphaFill: enriching AlphaFold models with ligand cofactors and ions. *Nat Methods* 2023;20:205–213. doi:10.1038/s41592-022-01685-y
- **GitHub**: https://github.com/PDB-REDO/alphafill
- **Applied**: cartoon + ligand stick + contact residue label layout; transparent pocket surface

## Implementation note
- PyMOL not assumed installed on the host system
- matplotlib schematic produced for publication figure (fig05_interaction.pdf/svg/png)
- Companion PyMOL script (fig05_interaction_pymol.pml) provided for users with PyMOL ≥ 2.5

## Data source
- Tier 1 validation table: `06_Report/Mr_Repro/results/benchmark/validation_table.tsv`
- Case: MCP-AlaIle, Vina = −5.61, CNN pose = 0.8995, CNN affinity = 4.416

## Style decisions
- Receptor: rainbow N→C (blue→red) cartoon helices
- Ligand: green (#00AA44) stick representation
- Pocket: sky blue ellipse, 18% alpha
- Contact residues: gray (#888888) sticks + labels
- Score annotation box at bottom
