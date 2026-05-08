"""
dad.core — Core algorithm modules for the DAD pipeline.

Module layout
-------------
structure   Stage 4 — ColabFold / AlphaFold2 structure prediction
pocket      Stage 5 — P2Rank pocket detection
docking     Stage 6-8 — ligand preparation, box config, GNINA docking
interaction Stage 9-10 — contact analysis, PLIP, aggregation & ranking
visualize   Stage 11 — py3Dmol HTML, ChimeraX cxc, report generation
triage      Stage 3 — pre-docking biological triage (Phobius/DeepTMHMM + HMMER)

All public functions are documented with NumPy-style docstrings.  The main
pipeline modules now contain executable Tier 1 replay and core workflow logic;
external-dataset ingest stubs live separately under Mr_Repro Phase D scope.
"""
