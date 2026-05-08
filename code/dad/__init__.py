"""
dad — Dynamic Affinity Dock
===========================

End-to-end protein-ligand many-to-many docking pipeline.

Stages
------
0  Project scope & positioning
1  Input ingestion (FASTA + SMILES)
2  Sequence QC & dereplication
3  Pre-docking biological triage  (dad.core.triage)
4  Structure prediction            (dad.core.structure)
5  Pocket detection                (dad.core.pocket)
6  Ligand preparation              (dad.core.docking)
7  Auto box configuration          (dad.core.docking)
8  GNINA docking                   (dad.core.docking)
9  Interaction profiling           (dad.core.interaction)
10 Aggregation & ranking           (dad.core.interaction)
11 Visualization & report          (dad.core.visualize)
12 Reproducibility manifest        (Snakemake + dad.core)

Target journal: Nature Protocols
"""

__version__ = "0.1.0"
__author__ = "DAD Team"
