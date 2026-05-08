"""
Figure 2 — Pre-docking biological triage decision tree
Reference style: PoseBusters Fig 1 (Buttenschoen et al. 2024 Chem Sci)
                 decision-tree node/branch layout
Input: 06_Report/Mr_Bio/rationale.md — 5 rules R1-R5
       Tier 1 case examples: NA23_RS01195=PASS_CLIPPED / NA23_RS08105=PASS
                             / NA23_RS00870=PASS_CLIPPED
Spec: 89 mm single-column, Arial, vertical tree
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
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

# colorblind-safe
COL = {
    "rule":    "#0072B2",
    "pass":    "#009E73",
    "clipped": "#56B4E9",
    "exclude": "#D55E00",
    "input":   "#555555",
    "example": "#CC79A7",
    "diamond": "#F0E442",
}

FIG_W_IN = 3.5   # 89 mm
FIG_H_IN = 7.2
fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
ax.set_xlim(0, 3.5)
ax.set_ylim(-0.3, 10.5)
ax.axis("off")

# ── Helper functions ──────────────────────────────────────────────────────────
def box(ax, cx, cy, w, h, text, fcolor, ecolor="white", tcolor="white",
        fontsize=6, bold=False, lw=0.8):
    rect = FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.04",
        facecolor=fcolor, edgecolor=ecolor, linewidth=lw, zorder=3,
        alpha=0.9,
    )
    ax.add_patch(rect)
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fontsize,
            color=tcolor, fontweight="bold" if bold else "normal",
            multialignment="center", zorder=4)

def diamond(ax, cx, cy, w, h, text, fcolor="#F0E442", tcolor="#222222",
            fontsize=6):
    pts = np.array([[cx, cy + h / 2], [cx + w / 2, cy],
                    [cx, cy - h / 2], [cx - w / 2, cy]])
    patch = plt.Polygon(pts, closed=True, facecolor=fcolor, edgecolor="white",
                        linewidth=0.8, zorder=3, alpha=0.92)
    ax.add_patch(patch)
    ax.text(cx, cy, text, ha="center", va="center", fontsize=fontsize,
            color=tcolor, fontweight="bold", multialignment="center", zorder=4)

def arrow(ax, x1, y1, x2, y2, label="", label_side="right",
          color="#555555"):
    ax.annotate(
        "",
        xy=(x2, y2),
        xytext=(x1, y1),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.7,
                        mutation_scale=7),
        zorder=2,
    )
    if label:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        dx = 0.1 if label_side == "right" else -0.1
        ax.text(mx + dx, my, label, ha="left" if label_side == "right"
                else "right", va="center", fontsize=5, color=color,
                fontstyle="italic")

def side_arrow(ax, x_from, y_from, x_to, y_to, label="", color="#555555"):
    ax.annotate(
        "",
        xy=(x_to, y_to),
        xytext=(x_from, y_from),
        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.7,
                        mutation_scale=7, connectionstyle="arc3,rad=0"),
        zorder=2,
    )
    if label:
        mx = (x_from + x_to) / 2
        my = (y_from + y_to) / 2
        ax.text(mx, my + 0.08, label, ha="center", va="bottom", fontsize=5,
                color=color, fontstyle="italic")

# ── Input box ────────────────────────────────────────────────────────────────
box(ax, 1.75, 10.0, 2.2, 0.42,
    "ORF Translation (FASTA input pool)",
    COL["input"], tcolor="white", fontsize=6, bold=True)

arrow(ax, 1.75, 9.78, 1.75, 9.42)

# ── R1 Length filter ─────────────────────────────────────────────────────────
diamond(ax, 1.75, 9.2, 2.4, 0.50,
        "R1  Length ≥ 50 aa?", fontsize=5.8)
arrow(ax, 1.75, 8.95, 1.75, 8.58, label="YES", color=COL["pass"])

# EXCLUDE: too short
side_arrow(ax, 2.95, 9.2, 3.3, 9.2, color=COL["exclude"])
box(ax, 3.3, 9.2, 0.35, 0.36, "EXCL\n(short)", COL["exclude"],
    fontsize=5)

# ── R2 Signal peptide ─────────────────────────────────────────────────────────
diamond(ax, 1.75, 8.35, 2.4, 0.50,
        "R2  Signal peptide\ndetected?", fontsize=5.8)
arrow(ax, 1.75, 8.1, 1.75, 7.73, label="NO", color="#555555")

# SP-only → PASS_CLIPPED branch
side_arrow(ax, 0.55, 8.35, 0.1, 8.35, color=COL["clipped"])
box(ax, 0.1, 7.78, 0.75, 0.85,
    "PASS\nCLIPPED\n(mature\nform)", COL["clipped"],
    tcolor="white", fontsize=5)
ax.text(0.46, 8.38, "YES →", ha="right", va="center", fontsize=5,
        color=COL["clipped"], fontstyle="italic")

# ── R3 nTM count ──────────────────────────────────────────────────────────────
diamond(ax, 1.75, 7.5, 2.4, 0.50,
        "R3  nTM helices?", fontsize=5.8)

# nTM=0 → PASS
arrow(ax, 1.75, 7.25, 1.75, 6.88, label="0", color=COL["pass"])
box(ax, 1.75, 6.65, 1.6, 0.42, "B: PASS (soluble)",
    COL["pass"], fontsize=5.5)
arrow(ax, 1.75, 6.44, 1.75, 6.08)

# nTM 1-2 branch (right)
arrow(ax, 2.95, 7.5, 3.2, 7.5, color=COL["clipped"])
box(ax, 3.2, 7.5, 0.55, 0.42,
    "C: CLIP\nextracell.", COL["clipped"], fontsize=5)
ax.text(2.96, 7.62, "1–2", ha="left", va="bottom", fontsize=5,
        color=COL["clipped"], fontstyle="italic")

# nTM 3-6 branch (left)
arrow(ax, 0.55, 7.5, 0.18, 7.5, color="#CC79A7")
box(ax, 0.1, 7.5, 0.54, 0.42,
    "D: CLIP\nloop≥60aa", "#CC79A7", fontsize=5)
ax.text(0.55, 7.62, "3–6", ha="right", va="bottom", fontsize=5,
        color="#CC79A7", fontstyle="italic")

# nTM ≥7 EXCLUDE
ax.text(1.75, 7.73, "≥7 → EXCLUDE (polytopic)", ha="center", va="bottom",
        fontsize=4.8, color=COL["exclude"], fontstyle="italic")
exc_box = FancyBboxPatch((0.82, 7.85), 1.86, 0.22,
                          boxstyle="round,pad=0.02",
                          facecolor=COL["exclude"], edgecolor="none",
                          alpha=0.15, zorder=1)
ax.add_patch(exc_box)

# ── R4 Dock-region ≥50 aa ─────────────────────────────────────────────────────
diamond(ax, 1.75, 5.85, 2.4, 0.50,
        "R4  Dock region ≥ 50 aa?", fontsize=5.8)
arrow(ax, 1.75, 5.6, 1.75, 5.22, label="YES", color=COL["pass"])
side_arrow(ax, 2.95, 5.85, 3.3, 5.85, color=COL["exclude"])
box(ax, 3.3, 5.85, 0.35, 0.36, "EXCL\n(short\npocket)", COL["exclude"],
    fontsize=5)
ax.text(2.95, 5.97, "NO", ha="left", va="bottom", fontsize=5,
        color=COL["exclude"], fontstyle="italic")

# ── R5 Functional class ───────────────────────────────────────────────────────
diamond(ax, 1.75, 5.0, 2.4, 0.50,
        "R5  Functional class\n(HMMER/Pfam)", fontsize=5.8)
arrow(ax, 1.75, 4.75, 1.75, 4.38)

# Boost
side_arrow(ax, 2.95, 5.0, 3.3, 5.0, color=COL["pass"])
box(ax, 3.3, 5.0, 0.35, 0.46,
    "+1 tier\nBOOST\n(MCP\nSBP)", COL["pass"], fontsize=5)
ax.text(2.96, 5.12, "sensor/SBP", ha="left", va="bottom", fontsize=4.8,
        color=COL["pass"], fontstyle="italic")

# Downrank DUF
side_arrow(ax, 0.55, 5.0, 0.18, 5.0, color="#888888")
box(ax, 0.1, 5.0, 0.35, 0.46,
    "FLAG\nDOWN-\nRANK\n(DUF)", "#888888", fontsize=5)
ax.text(0.54, 5.12, "DUF/hyp.", ha="right", va="bottom", fontsize=4.8,
        color="#888888", fontstyle="italic")

# ── Final output ──────────────────────────────────────────────────────────────
box(ax, 1.75, 4.15, 2.6, 0.42,
    "triage_report.tsv  →  Stage 4 (ColabFold)",
    COL["pass"], tcolor="white", fontsize=5.5, bold=True)

# ── Tier 1 case examples ──────────────────────────────────────────────────────
ax.text(1.75, 3.55, "Tier 1 examples", ha="center", va="top",
        fontsize=6, fontweight="bold", color=COL["example"])

examples = [
    ("NA23_RS01195  (MCP)", "nTM=2, SP−", "PASS_CLIPPED\n(periplasmic domain)"),
    ("NA23_RS08105  (CRP/FNR)", "nTM=0, SP−", "PASS\n(full soluble)"),
    ("NA23_RS00870  (RbsB)", "nTM=0, SP+", "PASS_CLIPPED\n(mature form)"),
]
y_ex = 3.35
for name, detail, outcome in examples:
    box(ax, 1.75, y_ex, 3.0, 0.46,
        f"{name}\n{detail}  →  {outcome}",
        COL["example"], tcolor="white", fontsize=5, lw=0.6)
    y_ex -= 0.56

# ── Panel label ───────────────────────────────────────────────────────────────
ax.text(0.02, 10.45, "b", fontsize=8, fontweight="bold", va="top", ha="left",
        transform=ax.transData)

fig.tight_layout(pad=0.3)

out_base = str(FIG_DIR / "fig02_triage_tree")
fig.savefig(out_base + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(out_base + ".svg", format="svg", bbox_inches="tight")
fig.savefig(out_base + ".png", format="png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 2 saved.")
