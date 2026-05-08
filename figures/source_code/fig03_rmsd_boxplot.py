"""
Figure 3 — RMSD distribution per protein family (16 RCSB seed redocking)
Reference style: PoseBusters Fig 4 grouped boxplot (Buttenschoen et al. 2024)
                 + DiffDock RMSD percentile plot (Corso et al. 2023)
Input: redocking_results.tsv (16 rows)
       redocking_summary.json (family aggregate)
Spec: 183mm double-column, 2-panel (top-pose vs best-of-9) grouped by family
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import csv, json
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
FIG_DIR    = SCRIPT_DIR.parent
DATA_DIR   = FIG_DIR.parent.parent / "06_Report" / "Mr_Repro" / \
             "external_validation" / "rcsb_seed"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ── Load data ─────────────────────────────────────────────────────────────────
DATA_PATH = DATA_DIR / "redocking_results.tsv"

rows = []
with open(DATA_PATH, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rows.append(row)

# Family display names
FAMILY_MAP = {
    "crp_fnr_regulator":       "CRP/FNR\nRegulator",
    "mcp_ligand_binding_domain": "MCP\nLBD",
    "amino_acid_binding":      "Amino-acid\nBinding",
    "periplasmic_sugar_binding": "Periplasmic\nSugar Binding",
}
FAMILY_ORDER = [
    "crp_fnr_regulator",
    "mcp_ligand_binding_domain",
    "amino_acid_binding",
    "periplasmic_sugar_binding",
]

# Group by family
family_data = {f: {"top": [], "best": [], "cases": []} for f in FAMILY_ORDER}
for r in rows:
    fam = r["family"]
    if fam in family_data:
        try:
            top_rmsd = float(r["top_pose_rmsd_a"])
            best_rmsd = float(r["best_pose_rmsd_a"])
        except ValueError:
            continue
        family_data[fam]["top"].append(top_rmsd)
        family_data[fam]["best"].append(best_rmsd)
        family_data[fam]["cases"].append({
            "case_id": r["case_id"],
            "top_rmsd": top_rmsd,
            "best_rmsd": best_rmsd,
            "redock_pass": r["redock_pass"],
            "best_pass": r["best_pose_pass"],
        })

# ── Colours ──────────────────────────────────────────────────────────────────
COL_TOP  = "#0072B2"   # blue — top pose
COL_BEST = "#56B4E9"   # sky blue — best-of-9 diagnostic
COL_FAIL = "#D55E00"   # vermillion — failure annotation
THRESH_COLOR = "#CC0000"

# ── Figure ───────────────────────────────────────────────────────────────────
FIG_W_IN = 7.2
FIG_H_IN = 3.6
fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

n_fam = len(FAMILY_ORDER)
group_w = 0.7
bar_w = 0.28
gap = 0.12
x_centers = np.arange(n_fam)

bp_top_list = []
bp_best_list = []

for i, fam in enumerate(FAMILY_ORDER):
    top_vals  = family_data[fam]["top"]
    best_vals = family_data[fam]["best"]

    x_top  = x_centers[i] - (bar_w / 2 + gap / 2)
    x_best = x_centers[i] + (bar_w / 2 + gap / 2)

    # jitter
    np.random.seed(42)
    jitter_top  = np.random.uniform(-0.06, 0.06, len(top_vals))
    jitter_best = np.random.uniform(-0.06, 0.06, len(best_vals))

    # box plot manually using matplotlib boxplot
    bp_top = ax.boxplot(
        top_vals,
        positions=[x_top],
        widths=bar_w,
        patch_artist=True,
        medianprops=dict(color="white", linewidth=1.2),
        boxprops=dict(facecolor=COL_TOP, alpha=0.85, linewidth=0.7),
        whiskerprops=dict(color=COL_TOP, linewidth=0.8),
        capprops=dict(color=COL_TOP, linewidth=0.8),
        flierprops=dict(marker="o", markersize=2.5, markerfacecolor=COL_TOP,
                        markeredgewidth=0.4, alpha=0.7),
        zorder=3,
    )
    bp_best = ax.boxplot(
        best_vals,
        positions=[x_best],
        widths=bar_w,
        patch_artist=True,
        medianprops=dict(color="white", linewidth=1.2),
        boxprops=dict(facecolor=COL_BEST, alpha=0.85, linewidth=0.7),
        whiskerprops=dict(color=COL_BEST, linewidth=0.8),
        capprops=dict(color=COL_BEST, linewidth=0.8),
        flierprops=dict(marker="o", markersize=2.5, markerfacecolor=COL_BEST,
                        markeredgewidth=0.4, alpha=0.7),
        zorder=3,
    )

    # Jittered dots
    ax.scatter([x_top + j for j in jitter_top], top_vals,
               color=COL_TOP, s=8, alpha=0.7, zorder=4, linewidths=0)
    ax.scatter([x_best + j for j in jitter_best], best_vals,
               color=COL_BEST, s=8, alpha=0.7, zorder=4, linewidths=0)

    # Annotate failure cases (top_rmsd > 2.0)
    for c in family_data[fam]["cases"]:
        if c["top_rmsd"] > 2.0:
            jj = np.random.uniform(-0.05, 0.05)
            ax.annotate(
                "",
                xy=(x_top + jj, c["top_rmsd"]),
                xytext=(x_top + jj, c["top_rmsd"] + 0.3),
                arrowprops=dict(arrowstyle="-|>", color=COL_FAIL, lw=0.7,
                                mutation_scale=5),
                zorder=5,
            )
            ax.scatter([x_top + jj], [c["top_rmsd"]],
                       color=COL_FAIL, s=14, zorder=6, linewidths=0.5,
                       edgecolors="white")

# ── 2.0 Å threshold line ──────────────────────────────────────────────────────
ax.axhline(2.0, color=THRESH_COLOR, linewidth=0.9, linestyle="--",
           zorder=2, label="2.0 Å threshold")
ax.text(n_fam - 0.1, 2.08, "2.0 Å", ha="right", va="bottom",
        fontsize=5.5, color=THRESH_COLOR)

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xticks(x_centers)
ax.set_xticklabels([FAMILY_MAP[f] for f in FAMILY_ORDER],
                   fontsize=6.5, multialignment="center")
ax.set_ylabel("RMSD to crystal pose (Å)", fontsize=7)
ax.set_xlabel("Protein family", fontsize=7)
ax.set_xlim(-0.6, n_fam - 0.4)
ax.set_ylim(-0.1, 14.5)
ax.set_title("RCSB Seed Redocking — RMSD by Protein Family\n"
             "(top-pose vs. diagnostic best-of-9)", fontsize=7, pad=4)

# ── Legend ────────────────────────────────────────────────────────────────────
leg_patches = [
    mpatches.Patch(facecolor=COL_TOP,  label="Top-pose RMSD (operational)",
                   alpha=0.85),
    mpatches.Patch(facecolor=COL_BEST, label="Best-of-9 RMSD (diagnostic)",
                   alpha=0.85),
    mpatches.Patch(facecolor=COL_FAIL, label="Top-pose FAIL (> 2.0 Å)",
                   alpha=0.85),
    mpatches.Patch(facecolor="none",   label="— 2.0 Å threshold",
                   linewidth=0.9, linestyle="--", edgecolor=THRESH_COLOR),
]
ax.legend(handles=leg_patches, fontsize=5.5, frameon=False,
          loc="upper left", bbox_to_anchor=(0.01, 0.99))

# ── n= annotation per family ─────────────────────────────────────────────────
for i, fam in enumerate(FAMILY_ORDER):
    n = len(family_data[fam]["top"])
    ax.text(x_centers[i], -0.7, f"n={n}", ha="center", va="top",
            fontsize=5.5, color="#555555")

# ── Panel label ───────────────────────────────────────────────────────────────
ax.text(-0.06, 1.02, "c", fontsize=8, fontweight="bold",
        transform=ax.transAxes, va="top", ha="left")

fig.tight_layout(pad=0.4)

out_base = str(FIG_DIR / "fig03_rmsd_boxplot")
fig.savefig(out_base + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(out_base + ".svg", format="svg", bbox_inches="tight")
fig.savefig(out_base + ".png", format="png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 3 saved.")
