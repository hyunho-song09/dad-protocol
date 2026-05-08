"""
Figure 5 — Interaction view: MCP-AlaIle best worked example
Reference style: ColabFold protocol Fig + AlphaFill protocol Fig
                 cartoon + ligand stick + binding residue label
Input: AW1_ref/structure_MCP.pdb (receptor cartoon, rainbow N->C)
       Tier 1 replay: Vina -5.61, CNN pose 0.8995
Approach: py3Dmol HTML screenshot / matplotlib schematic fallback
          (PyMOL not assumed installed)
Spec: 89mm single-column, receptor cartoon rainbow, ligand green stick,
      5Å contacts gray sticks + label, pocket surface transparent
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
FIG_DIR    = SCRIPT_DIR.parent

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 6.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

# ── NOTE: This figure uses a schematic representation ─────────────────────────
# PyMOL/ChimeraX structural render requires local installation.
# A companion PyMOL script (fig05_interaction_pymol.pml) is provided in
# source_code/ for users with PyMOL installed to generate the actual
# structural image. This matplotlib version serves as a publication-ready
# schematic figure with quantitative annotations.

# ── Color definitions ─────────────────────────────────────────────────────────
# Rainbow N->C coloring for receptor cartoon
rainbow_colors = ["#0000FF", "#0055FF", "#00AAFF", "#00FFAA",
                  "#55FF00", "#AAFF00", "#FFAA00", "#FF5500", "#FF0000"]
LIGAND_COLOR = "#00AA44"    # green for Ala-Ile stick
CONTACT_COLOR = "#888888"   # gray for contact residues
SURFACE_COLOR = "#56B4E9"   # sky blue pocket surface
POCKET_ALPHA  = 0.18

FIG_W_IN = 3.5
FIG_H_IN = 4.0
fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
ax.set_xlim(0, 3.5)
ax.set_ylim(0, 4.0)
ax.axis("off")
fig.patch.set_facecolor("#F8F8F8")
ax.set_facecolor("#F8F8F8")

# ── Background label ──────────────────────────────────────────────────────────
ax.text(1.75, 3.88, "MCP Ligand-Binding Domain — Ala-Ile Docking",
        ha="center", va="top", fontsize=6.5, fontweight="bold",
        color="#222222")
ax.text(1.75, 3.70, "Tier 1 worked example  |  Vina = −5.61 kcal/mol  "
        "|  CNN pose = 0.900",
        ha="center", va="top", fontsize=5.5, color="#444444")

# ── Cartoon representation of receptor helix bundle ───────────────────────────
# Draw stylized helices as thick curved ribbons (schematic)
np.random.seed(7)

# Periplasmic sensory domain schematic: 4 helices surrounding pocket
helix_centers = [
    (0.75, 2.3), (1.35, 2.75), (2.4, 2.75), (2.9, 2.3),
    (0.85, 1.6), (2.65, 1.6),
]
helix_lengths = [1.2, 1.0, 1.0, 1.2, 0.9, 0.9]
helix_angles  = [80, 100, 80, 100, 70, 110]
helix_colors  = [rainbow_colors[i % len(rainbow_colors)]
                 for i in range(len(helix_centers))]

for (hx, hy), hl, ha_deg, hcol in zip(helix_centers, helix_lengths,
                                        helix_angles, helix_colors):
    ha_rad = np.deg2rad(ha_deg)
    dx = hl * 0.5 * np.cos(ha_rad)
    dy = hl * 0.5 * np.sin(ha_rad)
    # Ribbon: thick curved tube
    ax.annotate(
        "",
        xy=(hx + dx, hy + dy),
        xytext=(hx - dx, hy - dy),
        arrowprops=dict(
            arrowstyle="-|>",
            color=hcol,
            lw=5.0,
            mutation_scale=12,
            alpha=0.82,
        ),
        zorder=3,
    )

# ── N- and C-terminus labels ──────────────────────────────────────────────────
ax.text(0.35, 3.12, "N", ha="center", va="center", fontsize=7,
        fontweight="bold", color=rainbow_colors[0],
        bbox=dict(boxstyle="circle,pad=0.15", facecolor=rainbow_colors[0],
                  edgecolor="white", linewidth=0.5, alpha=0.7))
ax.text(3.15, 1.25, "C", ha="center", va="center", fontsize=7,
        fontweight="bold", color=rainbow_colors[-1],
        bbox=dict(boxstyle="circle,pad=0.15", facecolor=rainbow_colors[-1],
                  edgecolor="white", linewidth=0.5, alpha=0.7))

# ── Pocket surface (transparent ellipse) ─────────────────────────────────────
pocket_ellipse = mpatches.Ellipse(
    (1.82, 2.18), 1.5, 1.2,
    angle=15,
    facecolor=SURFACE_COLOR, edgecolor=SURFACE_COLOR,
    alpha=POCKET_ALPHA, zorder=2,
)
ax.add_patch(pocket_ellipse)
ax.text(1.82, 2.72, "Binding pocket", ha="center", va="bottom",
        fontsize=5, color=SURFACE_COLOR, alpha=0.85, fontstyle="italic")

# ── Ligand Ala-Ile stick representation ──────────────────────────────────────
lig_cx, lig_cy = 1.82, 2.08
# Dipeptide backbone schematic: N-Cα-C=O-N-Cα-C
lig_atoms = [(lig_cx - 0.30, lig_cy + 0.05),   # N (Ala)
             (lig_cx - 0.10, lig_cy),            # Cα (Ala)
             (lig_cx + 0.12, lig_cy + 0.08),     # C=O
             (lig_cx + 0.12, lig_cy - 0.12),     # O
             (lig_cx + 0.28, lig_cy + 0.08),     # N (Ile)
             (lig_cx + 0.48, lig_cy),            # Cα (Ile)
             (lig_cx + 0.68, lig_cy + 0.08),     # Cβ (Ile)
             (lig_cx + 0.68, lig_cy - 0.15),     # Cγ
             ]

for j in range(len(lig_atoms) - 1):
    x0, y0 = lig_atoms[j]
    x1, y1 = lig_atoms[j + 1]
    lw = 2.8 if j != 3 else 1.8
    ax.plot([x0, x1], [y0, y1], color=LIGAND_COLOR, linewidth=lw,
            solid_capstyle="round", zorder=5)

for ax_pos, ay_pos in lig_atoms:
    ax.scatter([ax_pos], [ay_pos], color=LIGAND_COLOR, s=14,
               zorder=6, edgecolors="white", linewidths=0.4)

ax.text(lig_cx + 0.2, lig_cy - 0.30, "Ala-Ile",
        ha="center", va="top", fontsize=6, color=LIGAND_COLOR,
        fontweight="bold")

# ── 5 Å contact residues (gray sticks + labels) ──────────────────────────────
contact_residues = [
    ("Arg64",  1.20, 2.45),
    ("Thr67",  1.35, 1.78),
    ("Asn71",  1.55, 2.52),
    ("Trp130", 2.28, 2.48),
    ("Leu133", 2.55, 1.80),
    ("Ser137", 2.30, 1.65),
]

for res_name, rx, ry in contact_residues:
    # Stick from residue backbone to near pocket
    ax.plot([rx, lig_cx + np.random.uniform(-0.2, 0.2)],
            [ry, lig_cy + np.random.uniform(-0.15, 0.15)],
            color=CONTACT_COLOR, linewidth=1.0, linestyle="--",
            alpha=0.55, zorder=4)
    ax.scatter([rx], [ry], color=CONTACT_COLOR, s=10, zorder=5,
               edgecolors="white", linewidths=0.3)
    ax.text(rx, ry + 0.10, res_name, ha="center", va="bottom",
            fontsize=4.8, color=CONTACT_COLOR)

# ── Score annotation box ──────────────────────────────────────────────────────
score_box = FancyBboxPatch(
    (0.05, 0.05), 3.4, 0.55,
    boxstyle="round,pad=0.06",
    facecolor="white", edgecolor="#CCCCCC", linewidth=0.7, zorder=7,
)
ax.add_patch(score_box)
score_text = (
    "MCP (NA23_RS01195) × Ala-Ile\n"
    "Vina affinity: −5.61 kcal/mol  |  "
    "CNN pose: 0.900  |  CNN affinity: 4.416\n"
    "5 Å contacts: Arg64, Thr67, Asn71, Trp130, Leu133, Ser137"
)
ax.text(1.75, 0.33, score_text, ha="center", va="center",
        fontsize=5.0, color="#333333", multialignment="center", zorder=8)

# ── Colour bar (N→C rainbow) ─────────────────────────────────────────────────
cbar_x0, cbar_y0 = 0.1, 3.45
cbar_w, cbar_h = 1.8, 0.08
n_seg = len(rainbow_colors)
seg_w = cbar_w / n_seg
for seg_i, col in enumerate(rainbow_colors):
    rect = FancyBboxPatch(
        (cbar_x0 + seg_i * seg_w, cbar_y0), seg_w, cbar_h,
        boxstyle="square,pad=0",
        facecolor=col, edgecolor="none", zorder=4,
    )
    ax.add_patch(rect)
ax.text(cbar_x0, cbar_y0 + cbar_h + 0.04, "N", ha="center", va="bottom",
        fontsize=5.5, fontweight="bold", color=rainbow_colors[0])
ax.text(cbar_x0 + cbar_w, cbar_y0 + cbar_h + 0.04, "C", ha="center",
        va="bottom", fontsize=5.5, fontweight="bold", color=rainbow_colors[-1])
ax.text(cbar_x0 + cbar_w / 2, cbar_y0 + cbar_h + 0.04,
        "Receptor (rainbow N→C)", ha="center", va="bottom",
        fontsize=5, color="#444444")

# ── Legend ────────────────────────────────────────────────────────────────────
legend_items = [
    mpatches.Patch(facecolor=LIGAND_COLOR, label="Ala-Ile (ligand, sticks)"),
    mpatches.Patch(facecolor=CONTACT_COLOR, alpha=0.7,
                   label="5 Å contact residues"),
    mpatches.Patch(facecolor=SURFACE_COLOR, alpha=0.35,
                   label="Binding pocket (surface)"),
]
ax.legend(handles=legend_items, fontsize=4.8, frameon=True,
          framealpha=0.9, edgecolor="#DDDDDD",
          loc="lower right", bbox_to_anchor=(3.48, 0.64),
          bbox_transform=ax.transData)

# ── Note: schematic ───────────────────────────────────────────────────────────
ax.text(3.48, 1.48,
        "Schematic\n(PyMOL pml\navailable)",
        ha="right", va="top", fontsize=4.5, color="#999999",
        fontstyle="italic")

# ── Panel label ───────────────────────────────────────────────────────────────
ax.text(0.03, 3.98, "e", fontsize=8, fontweight="bold", va="top", ha="left",
        transform=ax.transData)

fig.tight_layout(pad=0.3)

out_base = str(FIG_DIR / "fig05_interaction")
fig.savefig(out_base + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(out_base + ".svg", format="svg", bbox_inches="tight")
fig.savefig(out_base + ".png", format="png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 5 saved.")
