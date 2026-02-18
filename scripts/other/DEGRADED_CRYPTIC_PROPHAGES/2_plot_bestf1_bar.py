#!/usr/bin/env python3
# Plot per-locus F1 scores (HQ vs LQ) as side-by-side bars across clustering levels in a 3x2 grid,
# with xtick labels for every K locus, sorted numerically (K1, KL2, K3, ...).
# NOW: columns = coverage (50, 80); rows = identity (00, 50, 80)

import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Iterable, Dict, Any, Optional, Set, List


def load_yaml_dict(path: str | Path) -> Dict[str, Any]:
    import yaml  # pip install pyyaml
    """Parse a YAML file and return its top-level mapping as a dict."""
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)  # safe loader avoids executing arbitrary tags
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise TypeError(f"YAML root is not a mapping; got {type(data).__name__}")
    return data


# ----------------------- USER CONFIG -----------------------
config_path = Path('../../../config/config.yaml')
config_dict = load_yaml_dict(config_path)

# input direcoties
user_path = config_dict['paths']['janusz']['main']
figshare_dir = Path(user_path, config_dict['paths']['janusz']['figshare_dir'])
supplement_dir = Path(user_path, config_dict['paths']['janusz']['supplement_dir'])

# output directories
FIGURES_DIR = Path().cwd() / "plots"
TABLES_DIR = Path().cwd() / "tables"

# output paths
in_path = Path(TABLES_DIR, "bestf1.tsv")
out_png = FIGURES_DIR / "bestf1_bar.png"
# out_pdf = FIGURES_DIR / "bestf1_bar.pdf"

Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)
Path(TABLES_DIR).mkdir(parents=True, exist_ok=True)

# Grid definition: columns = coverage (50, 80); rows = identity (00, 50, 80)
id_order  = [0, 50, 80]   # rows (top→bottom)
cov_order = [50, 80]      # columns (left→right)
required_versions = ["PCI00C50", "PCI00C80", "PCI50C50", "PCI50C80", "PCI80C50", "PCI80C80"]

# ---------- Load ----------
df = pd.read_csv(in_path, sep="\t", dtype=str)

# Required columns
for col in ["version", "locus", "quality"]:
    if col not in df.columns:
        raise ValueError(f"Missing required column: {col}")

# Numeric F1 helper (keep original if present)
if "F1_score_num" not in df.columns:
    if "F1_score" not in df.columns:
        raise ValueError("Need F1_score or F1_score_num column in the input.")
    df["F1_score_num"] = pd.to_numeric(df["F1_score"], errors="coerce")
else:
    df["F1_score_num"] = pd.to_numeric(df["F1_score_num"], errors="coerce")

# Keep only HQ/LQ
df = df[df["quality"].isin(["hq", "lq"])].copy()

# ---------- Build a global, numerically sorted K locus order ----------
# Parse e.g., "K1", "KL2", "K22", "KL102" -> 1, 2, 22, 102; non-matching go to the end.
def k_number_key(s: str):
    s = "" if s is np.nan else str(s)
    m = re.match(r"^\s*K(?:L)?\s*(\d+)\s*$", s, flags=re.IGNORECASE)
    return (int(m.group(1)) if m else float("inf"), s)

all_loci = sorted(df["locus"].dropna().unique(), key=k_number_key)

# ---------- Helpers ----------
def parse_ic(ver: str):
    # "PCI50C80" -> (identity=50, coverage=80)
    m = re.match(r"^PCI(\d{2})C(\d{2})$", str(ver))
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))

def grid_pos(ver: str):
    # SWAPPED: row by identity, column by coverage
    ident, cov = parse_ic(ver)
    if ident is None:
        return None
    try:
        r = id_order.index(ident)   # row = identity
        c = cov_order.index(cov)    # col = coverage
    except ValueError:
        return None
    return (r, c)

# ---------- Prepare figure ----------
fig, axes = plt.subplots(3, 2, figsize=(12, 6), sharex=True, sharey=True)

# Global y-limits: default to 0..1.05 (F1 is 0..1); widen if needed
all_vals = df["F1_score_num"].astype(float)
ymin = 0.0
ymax = max(1.05, float(np.nanmax(all_vals)) * 1.05 if np.isfinite(np.nanmax(all_vals)) else 1.05)

# ---------- Plot each panel ----------
bar_width = 0.4
x = np.arange(len(all_loci))  # shared x across panels

for ver in required_versions:
    pos = grid_pos(ver)
    if pos is None:
        continue
    r, c = pos
    ax = axes[r, c]

    sub = df[df["version"] == ver].copy()

    # Pivot to get columns 'hq' and 'lq' per locus and reindex to the global order
    pvt = (sub.pivot_table(index="locus", columns="quality", values="F1_score_num", aggfunc="max")
              .reindex(columns=["hq", "lq"]))
    pvt = pvt.reindex(all_loci)

    # Bars (missing values won't draw)
    ax.bar(x - bar_width/2, pvt["hq"].values, width=bar_width, label="HQ", align="center")
    ax.bar(x + bar_width/2, pvt["lq"].values, width=bar_width, label="LQ", align="center")

    ax.set_ylim(ymin, ymax)
    ax.axhline(0, lw=0.8, linestyle="-", alpha=0.4)

    # xticks for EVERY locus, sorted numerically
    ax.set_xticks(x)
    ax.set_xticklabels(all_loci, rotation=45, fontsize=7.5, fontweight='bold')
    ax.margins(x=0.005)



# Shared labels and headers
fig.supylabel("F1 score", x=0.05, fontsize=12, fontweight="bold")
fig.supxlabel("")
fig.suptitle("", y=0.98)

# # ---- Single shared labels + per-row/column values (identity on the RIGHT) ----
# # Big shared labels
# fig.text(0.5, 0.99, "Coverage", ha="center", va="top", fontsize=13, fontweight="bold")
# fig.text(0.995, 0.5, "Identity", ha="right", va="center",
#          rotation=-90, fontsize=13, fontweight="bold")

# Column values (above each coverage column)
for c2, cov in enumerate(cov_order):
    axes[0, c2].annotate(f"Coverage: {cov}%", xy=(0.5, 1.08), xycoords="axes fraction",
                         ha="center", va="bottom", fontsize=8, fontweight="bold")

# Row values (to the RIGHT of each identity row)
for r2, ident in enumerate(id_order):
    axes[r2, -1].annotate(f"Identity: {ident}%", xy=(1.02, 0.5), xycoords="axes fraction",
                          ha="left", va="center", fontsize=8, fontweight="bold")

for ax in axes.flat:
    plt.setp(ax.get_yticklabels(), fontweight="bold")


# One shared legend (top-right)
handles, labels = axes[0, 0].get_legend_handles_labels()
labels = ['High quality prophages', 'High and low quality prophages (default)']
fig.legend(handles, labels, loc="lower right", frameon=False)

fig.tight_layout(rect=[0.05, 0.05, 0.96, 0.92])
fig.savefig(out_png, dpi=300)
# fig.savefig(out_pdf)
print(f"Saved: {out_png}")
# print(f"Saved: {out_pdf}")
