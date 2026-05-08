# PyMOL script — DAD interaction view template (Figure 5 companion)
# Reference: ColabFold protocol fig + AlphaFill protocol fig style
# Run: pymol -c fig05_interaction_pymol.pml
#
# NOTE: This is a TEMPLATE script for generating three-dimensional structural
# renders using PyMOL (>= 2.5). It demonstrates the workflow using the RCSB
# co-crystal example TAR_ASP_2LIG as an EXAMPLE PDB ONLY — this is NOT a
# Tier 1 result and is NOT the MCP-AlaIle Tier 1 pair.
#
# To render the actual MCP-AlaIle Tier 1 interaction:
#   1. Produce a live docked SDF (Tier 1 is replay-only; no archived pose SDF).
#   2. Replace the receptor.pdb path below with the MCP structure.
#   3. Uncomment and update the ligand loading block.
#
# DAD_ROOT should be set to the project root, or edit paths below.

# ── Load structures (example: TAR_ASP_2LIG — EXAMPLE ONLY, not Tier 1) ───────
load D:/project/experiment/DAD/06_Report/Mr_Repro/external_validation/rcsb_seed/prepared/TAR_ASP_2LIG/receptor.pdb, receptor
# To load the docked ligand (replace path with actual docked SDF):
# load D:/project/experiment/DAD/06_Report/Mr_Repro/external_validation/rcsb_seed/redocking/TAR_ASP_2LIG/docked.sdf, ligand

# ── Display settings ──────────────────────────────────────────────────────────
bg_color white
set ray_opaque_background, 0
set antialias, 2
set ray_trace_mode, 1

# ── Receptor cartoon rainbow N->C ─────────────────────────────────────────────
hide everything, receptor
show cartoon, receptor
spectrum count, rainbow, receptor, byres=1
set cartoon_tube_radius, 0.25
set cartoon_fancy_helices, 1

# ── Ligand green sticks ───────────────────────────────────────────────────────
# show sticks, ligand
# color 0x00AA44, ligand
# set stick_radius, 0.15, ligand

# ── 5 Å contact residues ─────────────────────────────────────────────────────
# select contacts, byres (receptor within 5 of ligand)
# show sticks, contacts
# color gray70, contacts
# label contacts and name CA, "%s%s" % (resn, resi)

# ── Pocket surface ────────────────────────────────────────────────────────────
# create pocket_surf, byres (receptor within 8 of ligand)
# show surface, pocket_surf
# set transparency, 0.65, pocket_surf
# color 0x56B4E9, pocket_surf

# ── View and ray ─────────────────────────────────────────────────────────────
orient receptor
zoom receptor, 5
turn y, 20
turn x, -10

set ray_shadows, 0
set depth_cue, 0.5
ray 1800, 1800

png D:/project/experiment/DAD/Publication/figures/fig05_interaction_pymol_render.png, dpi=300
quit
