"""
Supplementary Figure S1 — Top-pose vs best-of-9 within-case ranking
Reference style: PoseBusters supplementary bar plot style
Input: redocking_results.tsv best_pose_index column (16 cases)
Codex R6 wording: "within-case ranking can place a bad pose first
                  even when good pose exists"
Spec: 183mm double-column, 16-case bar plot, x=case_id, y=best_pose_index
      color = redock_pass; index != 1 highlights ranking limitation
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import csv
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
    "xtick.labelsize": 5.8,
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

case_ids, best_indices, redock_pass_list, top_rmsds, best_rmsds = [], [], [], [], []

with open(DATA_PATH, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        try:
            case_ids.append(row["case_id"])
            best_indices.append(int(row["best_pose_index"]))
            redock_pass_list.append(row["redock_pass"])
            top_rmsds.append(float(row["top_pose_rmsd_a"]))
            best_rmsds.append(float(row["best_pose_rmsd_a"]))
        except (ValueError, KeyError):
            continue

n = len(case_ids)
x = np.arange(n)

# ── Colour by redock_pass AND by best_pose_index ──────────────────────────────
# Primary: PASS=green, FAIL=vermillion
# Secondary: if best_pose_index > 1 → add hatching to indicate ranking gap
COL_PASS   = "#009E73"
COL_FAIL   = "#D55E00"
COL_IDX1   = "#56B4E9"   # index=1, PASS (ideal)
HATCH_RANK = "///"       # ranking gap pattern

bar_colors = []
bar_hatches = []
for i in range(n):
    if redock_pass_list[i] == "PASS":
        bar_colors.append(COL_PASS)
    else:
        bar_colors.append(COL_FAIL)
    # Hatching if best pose is not index 1
    if best_indices[i] > 1:
        bar_hatches.append(HATCH_RANK)
    else:
        bar_hatches.append("")

# ── Figure ────────────────────────────────────────────────────────────────────
FIG_W_IN = 7.2
FIG_H_IN = 3.0
fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))

bars = ax.bar(x, best_indices, color=bar_colors, hatch=bar_hatches,
              edgecolor="white", linewidth=0.5, zorder=3, alpha=0.88,
              width=0.7)

# Horizontal reference at y=1 (index 1 = top pose = best pose)
ax.axhline(1, color="#888888", linewidth=0.8, linestyle="-", zorder=2)
ax.text(n - 0.5, 1.05, "index = 1\n(top pose IS best)", ha="right",
        va="bottom", fontsize=5.5, color="#888888", fontstyle="italic")

# Annotate top_pose_rmsd and best_pose_rmsd on each bar
for i in range(n):
    top_r = top_rmsds[i]
    best_r = best_rmsds[i]
    # Top RMSD label (inside bar top if space, else outside)
    bar_top = best_indices[i]
    label_y = bar_top + 0.15
    ax.text(x[i], label_y,
            f"top:{top_r:.2f}Å\nbest:{best_r:.2f}Å",
            ha="center", va="bottom", fontsize=4.0, color="#333333",
            multialignment="center")

# Best pose index = 1 annotation
index_ne_1 = [i for i in range(n) if best_indices[i] > 1]
for i in index_ne_1:
    ax.annotate(
        f"rank #{best_indices[i]}",
        xy=(x[i], best_indices[i]),
        xytext=(x[i], best_indices[i] + 1.2),
        ha="center", va="bottom", fontsize=4.8, color="#D55E00",
        arrowprops=dict(arrowstyle="-|>", color="#D55E00", lw=0.7,
                        mutation_scale=6),
        zorder=5,
    )

# ── Axes ──────────────────────────────────────────────────────────────────────
# Short case labels for readability
short_ids = []
for cid in case_ids:
    parts = cid.split("_")
    if len(parts) >= 2:
        short_ids.append(f"{parts[0]}\n{parts[1]}")
    else:
        short_ids.append(cid)

ax.set_xticks(x)
ax.set_xticklabels(short_ids, fontsize=5.2, rotation=0,
                   multialignment="center")
ax.set_yticks(range(1, 10))
ax.set_ylabel("Best-of-9 pose index (lowest RMSD)", fontsize=7)
ax.set_xlabel("Case (PDB / ligand)", fontsize=7)
ax.set_xlim(-0.6, n - 0.4)
ax.set_ylim(0, 11)
ax.set_title("Supplementary Figure S1 — Within-case pose ranking\n"
             "(best-of-9 index; index > 1 indicates CNN ranking gap)",
             fontsize=7, pad=4)

# ── Legend ────────────────────────────────────────────────────────────────────
leg_patches = [
    mpatches.Patch(facecolor=COL_PASS, label="Top-pose PASS (RMSD ≤ 2.0 Å)"),
    mpatches.Patch(facecolor=COL_FAIL, label="Top-pose FAIL (RMSD > 2.0 Å)"),
    mpatches.Patch(facecolor="#CCCCCC", hatch=HATCH_RANK,
                   label="Best-of-9 index > 1 (ranking gap — diagnostic)"),
]
ax.legend(handles=leg_patches, fontsize=5.5, frameon=False, loc="upper left")

# ── Value statement at bottom ─────────────────────────────────────────────────
n_rank_gap = len(index_ne_1)
ax.text(0.5, -0.18,
        f"{n_rank_gap}/{n} cases have best-of-9 index > 1, directly evidencing "
        f"CNN top-pose ranking limitation (diagnostic; not an operational metric).",
        ha="center", va="top", fontsize=5.2, color="#444444",
        fontstyle="italic", transform=ax.transAxes)

# ── Panel label ───────────────────────────────────────────────────────────────
ax.text(-0.06, 1.04, "S1", fontsize=7, fontweight="bold",
        transform=ax.transAxes, va="bottom", ha="left")

fig.tight_layout(pad=0.5)

out_base = str(FIG_DIR / "figS01_pose_ranking")
fig.savefig(out_base + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(out_base + ".svg", format="svg", bbox_inches="tight")
fig.savefig(out_base + ".png", format="png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Supplementary Figure S1 saved.")
