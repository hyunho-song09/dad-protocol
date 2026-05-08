"""
Figure 1 — DAD 12-stage Workflow Schematic
Reference style: ColabFold Nature Methods Fig 1 (Mirdita et al. 2022)
                 + AlphaFill Nature Methods Fig 1 (Hekkelman et al. 2022)
                 + Nature Protocols standard stage flow schematic
Spec: 183 mm double-column, Arial font, colorblind-safe phase palette
      Horizontal layout left->right, 4 phases color-coded
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
FIG_DIR    = SCRIPT_DIR.parent

# ── Nature journal standards ──────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,   # embed fonts in PDF
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

# ── Phase colour palette (colorblind-safe, adapted from Wong 2011 palette) ────
PHASE_COLORS = {
    "input":     "#0072B2",   # blue
    "triage":    "#D55E00",   # vermillion
    "structure": "#009E73",   # green
    "dock":      "#CC79A7",   # pink/magenta
    "output":    "#56B4E9",   # sky blue
    "decision":  "#F0E442",   # yellow (decision diamond)
}
PHASE_ALPHA = 0.85
BOX_TEXT_COLOR = "white"
DECISION_TEXT_COLOR = "#222222"

# ── Stage definitions ─────────────────────────────────────────────────────────
# Each tuple: (stage_num, label_line1, label_line2, phase)
STAGES = [
    (0,  "Input\nIngestion",           "Stage 1",  "input"),
    (1,  "Sequence QC &\nDereplicate", "Stage 2",  "input"),
    (2,  "Biological\nTriage",         "Stage 3",  "triage"),   # decision
    (3,  "Structure\nPrediction",      "Stage 4",  "structure"),
    (4,  "Pocket\nDetection",          "Stage 5",  "structure"),
    (5,  "Ligand\nPreparation",        "Stage 6",  "structure"),
    (6,  "Auto Box\nConfig",           "Stage 7",  "dock"),
    (7,  "Docking\n(GNINA 1.3)",       "Stage 8",  "dock"),
    (8,  "Interaction\nProfiling",     "Stage 9",  "dock"),
    (9,  "Aggregation\n& Ranking",     "Stage 10", "dock"),
    (10, "Visualization\n& Report",    "Stage 11", "output"),
]

PHASE_LABELS = [
    ("INPUT", "input"),
    ("TRIAGE", "triage"),
    ("STRUCTURE &\nLIGAND PREP", "structure"),
    ("DOCKING &\nANALYSIS", "dock"),
    ("OUTPUT", "output"),
]

# ── Canvas ────────────────────────────────────────────────────────────────────
FIG_W_IN = 7.2   # 183 mm in inches
FIG_H_IN = 3.0
fig, ax = plt.subplots(figsize=(FIG_W_IN, FIG_H_IN))
ax.set_xlim(0, 11.5)
ax.set_ylim(-0.3, 2.5)
ax.axis("off")

# ── Helper: rounded box ───────────────────────────────────────────────────────
BOX_W = 0.85
BOX_H = 0.60
BOX_Y = 0.9

def draw_stage_box(ax, cx, cy, label, stage_num, phase, is_decision=False):
    color = PHASE_COLORS[phase]
    if is_decision:
        # diamond
        diamond = mpatches.FancyArrowPatch(
            posA=(cx - 0.45, cy),
            posB=(cx + 0.45, cy),
            arrowstyle="-",
            color=color,
        )
        pts = np.array([[cx, cy + 0.33], [cx + 0.48, cy],
                        [cx, cy - 0.33], [cx - 0.48, cy]])
        diamond_patch = plt.Polygon(pts, closed=True, facecolor=color,
                                    edgecolor="white", linewidth=0.8,
                                    zorder=3, alpha=PHASE_ALPHA)
        ax.add_patch(diamond_patch)
        ax.text(cx, cy + 0.05, label, ha="center", va="center",
                fontsize=5.5, color=DECISION_TEXT_COLOR, fontweight="bold",
                zorder=4, multialignment="center")
        ax.text(cx, cy - 0.22, stage_num, ha="center", va="center",
                fontsize=4.5, color=DECISION_TEXT_COLOR, zorder=4)
    else:
        box = FancyBboxPatch(
            (cx - BOX_W / 2, cy - BOX_H / 2),
            BOX_W, BOX_H,
            boxstyle="round,pad=0.04",
            facecolor=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
            alpha=PHASE_ALPHA,
        )
        ax.add_patch(box)
        ax.text(cx, cy + 0.06, label, ha="center", va="center",
                fontsize=5.5, color=BOX_TEXT_COLOR, fontweight="bold",
                zorder=4, multialignment="center")
        ax.text(cx, cy - 0.19, stage_num, ha="center", va="center",
                fontsize=4.5, color=BOX_TEXT_COLOR, alpha=0.85, zorder=4)

def draw_arrow(ax, x1, x2, y, color="#555555"):
    ax.annotate(
        "",
        xy=(x2 - BOX_W / 2 + 0.01, y),
        xytext=(x1 + BOX_W / 2 - 0.01, y),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            lw=0.8,
            mutation_scale=8,
        ),
        zorder=2,
    )

# ── Phase background bands ────────────────────────────────────────────────────
phase_spans = [
    (0.0, 1.95,  "input"),
    (1.95, 2.45, "triage"),
    (2.45, 6.45, "structure"),
    (6.45, 10.45, "dock"),
    (10.45, 11.5, "output"),
]

for x0, x1, phase in phase_spans:
    band = FancyBboxPatch(
        (x0, BOX_Y - BOX_H / 2 - 0.12),
        x1 - x0,
        BOX_H + 0.24,
        boxstyle="round,pad=0.05",
        facecolor=PHASE_COLORS[phase],
        edgecolor="none",
        alpha=0.12,
        zorder=1,
    )
    ax.add_patch(band)

# ── User persona entry at far left ───────────────────────────────────────────
persona_x = -0.05
persona_texts = ["FASTA\n(ORF pool)", "SMILES\n(metabolites)"]
for i, txt in enumerate(persona_texts):
    yy = BOX_Y + 0.22 - i * 0.44
    ax.annotate(
        txt,
        xy=(0.08, BOX_Y),
        xytext=(persona_x, yy),
        fontsize=5.5,
        ha="right",
        va="center",
        color="#333333",
        arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.6,
                        mutation_scale=6),
        zorder=5,
    )

# ── Draw stages ───────────────────────────────────────────────────────────────
x_positions = np.linspace(0.5, 10.8, len(STAGES))

for i, (idx, label, stage_num, phase) in enumerate(STAGES):
    cx = x_positions[i]
    is_decision = (stage_num == "Stage 3")
    draw_stage_box(ax, cx, BOX_Y, label, stage_num, phase,
                   is_decision=is_decision)
    if i > 0:
        draw_arrow(ax, x_positions[i - 1], cx, BOX_Y,
                   color=PHASE_COLORS[phase])

# ── DeepTMHMM side branch from Stage 3 ──────────────────────────────────────
triage_x = x_positions[2]
exclude_y = BOX_Y - 0.75
exclude_box = FancyBboxPatch(
    (triage_x - 0.5, exclude_y - 0.18),
    1.0, 0.36,
    boxstyle="round,pad=0.03",
    facecolor="#E8E8E8",
    edgecolor="#AAAAAA",
    linewidth=0.6,
    zorder=3,
)
ax.add_patch(exclude_box)
ax.text(triage_x, exclude_y, "EXCLUDE\n(nTM ≥ 7)", ha="center",
        va="center", fontsize=5, color="#666666", zorder=4,
        multialignment="center")
ax.annotate(
    "",
    xy=(triage_x, exclude_y + 0.18),
    xytext=(triage_x, BOX_Y - 0.33),
    arrowprops=dict(arrowstyle="-|>", color="#888888", lw=0.7,
                    mutation_scale=7),
    zorder=2,
)
ax.text(triage_x + 0.26, BOX_Y - 0.55, "FAIL", ha="left", va="center",
        fontsize=5, color="#D55E00", fontstyle="italic")

# ── Output annotation at far right ───────────────────────────────────────────
out_x = x_positions[-1] + 0.55
ax.annotate(
    "N×M\nRanking\nMatrix",
    xy=(x_positions[-1] + BOX_W / 2, BOX_Y),
    xytext=(out_x, BOX_Y),
    fontsize=5.5,
    ha="left",
    va="center",
    color="#0072B2",
    fontweight="bold",
    arrowprops=dict(arrowstyle="-|>", color="#0072B2", lw=0.7,
                    mutation_scale=7),
    zorder=5,
)

# ── Phase label row (top) ─────────────────────────────────────────────────────
phase_label_configs = [
    (0.95, "INPUT", "input"),
    (2.2,  "TRIAGE", "triage"),
    (4.45, "STRUCTURE &\nLIGAND PREP", "structure"),
    (8.45, "DOCKING &\nANALYSIS", "dock"),
    (10.8, "OUTPUT", "output"),
]
for lx, ltxt, lphase in phase_label_configs:
    ax.text(lx, 1.62, ltxt, ha="center", va="bottom", fontsize=5.5,
            color=PHASE_COLORS[lphase], fontweight="bold",
            multialignment="center")

# ── PASS_CLIPPED annotation below Stage 3 ────────────────────────────────────
ax.text(triage_x, BOX_Y - 1.02,
        "PASS_CLIPPED: clip TM domain,\ndock periplasmic region",
        ha="center", va="top", fontsize=4.8, color="#555555",
        fontstyle="italic", multialignment="center")

# ── Panel label ───────────────────────────────────────────────────────────────
ax.text(-0.5, 2.45, "a", fontsize=8, fontweight="bold", va="top", ha="left",
        transform=ax.transData)

fig.tight_layout(pad=0.3)

out_base = str(FIG_DIR / "fig01_workflow")
fig.savefig(out_base + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(out_base + ".svg", format="svg", bbox_inches="tight")
fig.savefig(out_base + ".png", format="png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Figure 1 saved.")
