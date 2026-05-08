"""
Figure 4 — CNN score calibration + ROC curve
Reference style: GNINA paper score calibration scatter (McNutt et al. 2021)
                 standard ROC curve with AUROC annotation
Input: redocking_results.tsv (cnn_pose_score, top_pose_rmsd_a, redock_pass)
Spec: 183mm double-column, 2-panel (a) scatter (b) ROC

All metrics (Pearson r, AUROC, AUPRC, EF@25%) are computed directly from
the TSV at runtime — no values are hard-coded.
Computed values are written to:
  ../source_data/fig04_metrics.csv
  ../source_data/fig04_source_data.csv
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import csv
import os
from pathlib import Path
from datetime import datetime, timezone
from scipy.stats import pearsonr, spearmanr

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

# ── Paths (relative to this script's location) ────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
FIG_DIR    = SCRIPT_DIR.parent
DATA_DIR   = FIG_DIR.parent.parent / "06_Report" / "Mr_Repro" / \
             "external_validation" / "rcsb_seed"
DATA_PATH  = DATA_DIR / "redocking_results.tsv"
SOURCE_DATA_DIR = FIG_DIR / "source_data"
SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
case_ids, cnn_scores, top_rmsds, best_rmsds, redock_pass_flags = \
    [], [], [], [], []

with open(DATA_PATH, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        try:
            case_ids.append(row["case_id"])
            cnn_scores.append(float(row["cnn_pose_score"]))
            top_rmsds.append(float(row["top_pose_rmsd_a"]))
            best_rmsds.append(float(row["best_pose_rmsd_a"]))
            redock_pass_flags.append(1 if row["redock_pass"] == "PASS" else 0)
        except (ValueError, KeyError):
            continue

cnn_scores = np.array(cnn_scores)
top_rmsds  = np.array(top_rmsds)
best_rmsds = np.array(best_rmsds)
labels     = np.array(redock_pass_flags)   # 1=PASS, 0=FAIL
n = len(labels)

# ── Write source_data/fig04_source_data.csv ───────────────────────────────────
src_rows = []
for i in range(n):
    src_rows.append({
        "case_id":         case_ids[i],
        "cnn_pose_score":  f"{cnn_scores[i]:.10f}",
        "top_pose_rmsd_a": f"{top_rmsds[i]:.3f}",
        "best_pose_rmsd_a": f"{best_rmsds[i]:.3f}",
        "redock_pass":     "PASS" if labels[i] == 1 else "FAIL",
    })
src_path = SOURCE_DATA_DIR / "fig04_source_data.csv"
with open(src_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(src_rows[0].keys()))
    writer.writeheader()
    writer.writerows(src_rows)

# ── Compute Pearson r (top-pose vs CNN, best-of-9 vs CNN) ────────────────────
r_top,  p_top  = pearsonr(cnn_scores, top_rmsds)
r_best, p_best = pearsonr(cnn_scores, best_rmsds)
rho_top,  _    = spearmanr(cnn_scores, top_rmsds)
rho_best, _    = spearmanr(cnn_scores, best_rmsds)

# ── Empirical ROC ─────────────────────────────────────────────────────────────
order      = np.argsort(cnn_scores)[::-1]
lbl_sorted = labels[order]

n_pos = int(labels.sum())
n_neg = n - n_pos

tpr_list, fpr_list = [0.0], [0.0]
tp, fp = 0, 0
for lbl in lbl_sorted:
    if lbl == 1:
        tp += 1
    else:
        fp += 1
    tpr_list.append(tp / n_pos if n_pos else 0.0)
    fpr_list.append(fp / n_neg if n_neg else 0.0)
tpr_list.append(1.0)
fpr_list.append(1.0)

tpr_arr = np.array(tpr_list)
fpr_arr = np.array(fpr_list)
auroc   = float(np.trapezoid(tpr_arr, fpr_arr))

# ── Empirical AUPRC (Precision-Recall) ───────────────────────────────────────
prec_list, rec_list = [], []
tp_pr = 0
precision_sum = 0.0
for k, lbl in enumerate(lbl_sorted, start=1):
    if lbl == 1:
        tp_pr += 1
        precision_sum += tp_pr / k
    prec_list.append(tp_pr / k)
    rec_list.append(tp_pr / n_pos if n_pos else 0.0)

# Average precision matches 06_Report/Mr_Repro/benchmark/enrichment_metrics.py.
auprc = float(precision_sum / n_pos) if n_pos else 0.0

# ── Enrichment Factor at 25% ──────────────────────────────────────────────────
k25     = max(1, int(np.ceil(0.25 * n)))
top_lbl = lbl_sorted[:k25]
ef25    = float((top_lbl.sum() / k25) / (n_pos / n)) if n_pos else 0.0

# ── Write source_data/fig04_metrics.csv ──────────────────────────────────────
metrics = {
    "n":                         n,
    "n_pass_top_pose":           n_pos,
    "n_fail_top_pose":           n_neg,
    "pearson_r_top_vs_cnn":      round(r_top,  6),
    "pearson_p_top_vs_cnn":      round(p_top,  6),
    "pearson_r_best9_vs_cnn":    round(r_best, 6),
    "pearson_p_best9_vs_cnn":    round(p_best, 6),
    "spearman_rho_top_vs_cnn":   round(rho_top,  6),
    "spearman_rho_best9_vs_cnn": round(rho_best, 6),
    "auroc":                     round(auroc,  6),
    "auprc":                     round(auprc,  6),
    "ef_at_25pct":               round(ef25,   6),
    "rmsd_threshold_a":          2.0,
    "source_tsv":                str(DATA_PATH),
    "script":                    str(Path(__file__).resolve()),
    "generated_at":              datetime.now(timezone.utc)
                                 .strftime("%Y-%m-%dT%H:%M:%SZ"),
}
met_path = SOURCE_DATA_DIR / "fig04_metrics.csv"
with open(met_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["metric", "value"])
    for k, v in metrics.items():
        writer.writerow([k, v])

print(f"Computed metrics:")
print(f"  n={n}, n_pos={n_pos}, n_neg={n_neg}")
print(f"  Pearson r (top vs CNN) = {r_top:.6f}  p={p_top:.4g}")
print(f"  Pearson r (best9 vs CNN)= {r_best:.6f}  p={p_best:.4g}")
print(f"  AUROC  = {auroc:.6f}")
print(f"  AUPRC  = {auprc:.6f}")
print(f"  EF@25% = {ef25:.6f}")

# ── Colour palette ────────────────────────────────────────────────────────────
COL_PASS = "#009E73"
COL_FAIL = "#D55E00"
COL_ROC  = "#0072B2"
COL_REF  = "#BBBBBB"

# ── Figure ────────────────────────────────────────────────────────────────────
FIG_W_IN = 7.2
FIG_H_IN = 3.4
fig = plt.figure(figsize=(FIG_W_IN, FIG_H_IN))
gs  = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1], wspace=0.35)

ax_scatter = fig.add_subplot(gs[0])
ax_roc     = fig.add_subplot(gs[1])

# ── Panel (a): CNN pose score vs top-pose RMSD scatter ───────────────────────
pass_mask = labels == 1
fail_mask = labels == 0

ax_scatter.scatter(cnn_scores[pass_mask], top_rmsds[pass_mask],
                   color=COL_PASS, s=28, alpha=0.85, zorder=4,
                   label="PASS (top-pose RMSD ≤ 2.0 Å)",
                   edgecolors="white", linewidths=0.4)
ax_scatter.scatter(cnn_scores[fail_mask], top_rmsds[fail_mask],
                   color=COL_FAIL, s=28, alpha=0.85, zorder=4,
                   label="FAIL (top-pose RMSD > 2.0 Å)",
                   edgecolors="white", linewidths=0.4, marker="^")

# Regression line
m_fit, b_fit = np.polyfit(cnn_scores, top_rmsds, 1)
x_line = np.linspace(cnn_scores.min() - 0.01, cnn_scores.max() + 0.01, 100)
ax_scatter.plot(x_line, m_fit * x_line + b_fit, color="#555555", linewidth=1.0,
                linestyle="--", zorder=3)

# 2.0 Å threshold
ax_scatter.axhline(2.0, color="#CC0000", linewidth=0.8, linestyle=":",
                   zorder=2)
ax_scatter.text(0.88, 2.1, "2.0 Å", ha="left", va="bottom", fontsize=5.5,
                color="#CC0000", transform=ax_scatter.get_yaxis_transform())

# Pearson r annotation (computed from data)
ax_scatter.text(0.97, 0.97,
                f"Pearson r = {r_top:.3f}",
                transform=ax_scatter.transAxes,
                ha="right", va="top", fontsize=6,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor="#CCCCCC", linewidth=0.6))

ax_scatter.set_xlabel("CNN pose score", fontsize=7)
ax_scatter.set_ylabel("Top-pose RMSD to crystal pose (Å)", fontsize=7)
ax_scatter.set_title("(a) CNN score vs. top-pose RMSD\n"
                     f"(n = {n}; {n_pos} PASS, {n_neg} FAIL)",
                     fontsize=7, pad=4)
ax_scatter.legend(fontsize=5.5, frameon=False, loc="upper left")
ax_scatter.text(-0.12, 1.03, "a", fontsize=8, fontweight="bold",
                transform=ax_scatter.transAxes, va="bottom", ha="left")

# ── Panel (b): ROC curve ──────────────────────────────────────────────────────
ax_roc.plot(fpr_arr, tpr_arr, color=COL_ROC, linewidth=1.5,
            label=f"GNINA CNN score (AUROC = {auroc:.3f})", zorder=4)
ax_roc.plot([0, 1], [0, 1], color=COL_REF, linewidth=0.8,
            linestyle="--", zorder=2, label="Random (AUROC = 0.5)")

# Metrics text box — all computed from data
ax_roc.text(0.97, 0.08,
            f"AUROC = {auroc:.3f}\nAUPRC = {auprc:.3f}\nEF@25% = {ef25:.3f}",
            transform=ax_roc.transAxes,
            ha="right", va="bottom", fontsize=6,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#CCCCCC", linewidth=0.6))

ax_roc.set_xlabel("False Positive Rate", fontsize=7)
ax_roc.set_ylabel("True Positive Rate", fontsize=7)
ax_roc.set_title("(b) ROC — top-pose success\n"
                 "(cross-case classifier, n = 16 seed)",
                 fontsize=7, pad=4)
ax_roc.set_xlim(-0.02, 1.02)
ax_roc.set_ylim(-0.02, 1.05)
ax_roc.set_aspect("equal")
ax_roc.legend(fontsize=5.5, frameon=False, loc="lower right")
ax_roc.text(-0.15, 1.03, "b", fontsize=8, fontweight="bold",
            transform=ax_roc.transAxes, va="bottom", ha="left")

fig.tight_layout(pad=0.4)

out_base = str(FIG_DIR / "fig04_cnn_roc")
fig.savefig(out_base + ".pdf", format="pdf", bbox_inches="tight")
fig.savefig(out_base + ".svg", format="svg", bbox_inches="tight")
fig.savefig(out_base + ".png", format="png", dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Figure 4 saved to {out_base}.*")
print(f"Source data written to {SOURCE_DATA_DIR}/")
