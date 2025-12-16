#!/usr/bin/env python3
"""
Compact UpSet-style comparison of ST membership across three datasets (exclusive subsets).

- Top: black bars with counts (y-ticks every 500: "ST counts")
- Bottom-left: mini horizontal bars with per-dataset totals (bars grow leftward),
  labels positioned per-bar using data coords; Y order matches matrix (via invert_yaxis)
- Bottom-right: membership matrix (black=member, light gray=non-member; vertical line for intersections)
- Row order (matrix, top->bottom): KlebNNSsero, GWAS, Pathogenwatch
- Columns: exclusive subsets in custom order; empty subsets omitted

Input required columns: dataset, genomeID, ST, country
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def load_yaml_dict(path: str | Path) -> Dict[str, Any]:
    import yaml  # pip install pyyaml
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise TypeError(f"YAML root is not a mapping; got {type(data).__name__}")
    return data

# ----------------------- USER CONFIG -----------------------
CONFIG_PATH = Path('../../../config/config.yaml')
CONFIG_DICT = load_yaml_dict(CONFIG_PATH)

FIGSHARE_DIR = Path(CONFIG_DICT['paths']['janusz']['main']) / Path(CONFIG_DICT['paths']['janusz']['figshare_data']) / 'REVIEW'
OUTPUT_DIR   = Path(CONFIG_DICT['paths']['janusz']['main']) / Path(CONFIG_DICT['paths']['janusz']['output']) / "STS_GEODISTRIBUTION"
PLOTS_DIR    = OUTPUT_DIR / "plots"
TABLES_DIR   = OUTPUT_DIR / "tables"

# ----------------------- FILES -----------------------
IN_PATH = Path(TABLES_DIR, "metadata.tsv")
OUT_PNG = Path(PLOTS_DIR, "st_upset.png")

# ----------------------- DATASET LABELS -----------------------
A_LABEL = "KlebAMRnet"         # Pathogenwatch
B_LABEL = "KlebNNSsero"        # KlebNNSsero
C_LABEL = "klebpavia+kaspah"   # GWAS

DATASET_DISPLAY: Dict[str, str] = {
    A_LABEL: "KlebAMRnet",
    B_LABEL: "KlebNNSsero",
    C_LABEL: "GWAS",
}

# Matrix row order (top->bottom)
ROW_ORDER = [B_LABEL, C_LABEL, A_LABEL]  # KlebNNSsero, GWAS, Pathogenwatch

# Exclusive subsets and column order
SUBSET_ORDER = ["A unique", "A∩C only", "C unique", "A∩B only", "A∩B∩C", "B unique"]
SUBSET_MEMBERS = {
    "A unique":  {A_LABEL},
    "A∩C only":  {A_LABEL, C_LABEL},
    "C unique":  {C_LABEL},
    "A∩B only":  {A_LABEL, B_LABEL},
    "A∩B∩C":     {A_LABEL, B_LABEL, C_LABEL},
    "B unique":  {B_LABEL},
}

# ----------------------- STYLE (compact) -----------------------
FIGSIZE     = (5, 4)   # ↓ smaller, but fonts unchanged
BAR_COLOR   = "black"
BAR_EDGE    = "black"
BAR_LW      = 0.4
DOT_ON_FILL  = "black"
DOT_OFF_FILL = "#D9D9D9"
DOT_EDGE     = "black"
DOT_SIZE     = 180         # slightly smaller dots to reduce visual crowding
LINE_COLOR   = "black"
LINE_WIDTH   = 3
ROW_SHADE    = "#F0F0F0"

# ----------------------- IO PREP -----------------------
PLOTS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------- HELPERS -----------------------
def load_metadata(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    need = {"dataset", "genomeID", "ST", "country"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    return df

def build_sets(df: pd.DataFrame) -> Tuple[Set[str], Set[str], Set[str]]:
    S_A = set(df.loc[df["dataset"] == A_LABEL, "ST"].dropna().unique())
    S_B = set(df.loc[df["dataset"] == B_LABEL, "ST"].dropna().unique())
    S_C = set(df.loc[df["dataset"] == C_LABEL, "ST"].dropna().unique())
    present = set(df["dataset"].unique())
    expected = {A_LABEL, B_LABEL, C_LABEL}
    if not expected.issubset(present):
        raise ValueError(f"Expected datasets {sorted(expected)} not all present; observed {sorted(present)}")
    return S_A, S_B, S_C

def compute_exclusive_subsets(S_A: Set[str], S_B: Set[str], S_C: Set[str]) -> Dict[str, Set[str]]:
    subsets = {
        "A unique": (S_A - (S_B | S_C)),
        "A∩C only": ((S_A & S_C) - S_B),
        "C unique": (S_C - (S_A | S_B)),
        "A∩B only": ((S_A & S_B) - S_C),
        "A∩B∩C":    (S_A & S_B & S_C),
        "B unique": (S_B - (S_A | S_C)),
    }
    return {k: subsets[k] for k in SUBSET_ORDER if len(subsets[k]) > 0}

# ----------------------- PLOTTING -----------------------
def plot_upset(subsets: Dict[str, Set[str]], totals_by_ds: Dict[str, int]) -> None:
    x_keys = list(subsets.keys())
    n_cols = len(x_keys)
    counts = [len(subsets[k]) for k in x_keys]

    fig = plt.figure(figsize=FIGSIZE)
    gs = fig.add_gridspec(
        nrows=3, ncols=2,
        height_ratios=[2, 0.06, 1.25],   # ↓ tighter vertical allocation
        width_ratios=[0.9, 6.55],          # ↓ smaller left panel
        hspace=0.02, wspace=1            # ↓ tighter inter-axes spacing
    )
    ax_top  = fig.add_subplot(gs[0, 1])  # top bars
    ax_left = fig.add_subplot(gs[2, 0])  # left mini-bars
    ax_mat  = fig.add_subplot(gs[2, 1])  # membership matrix

    # ----- Top bars -----
    bars = ax_top.bar(range(n_cols), counts, color=BAR_COLOR, edgecolor=BAR_EDGE, linewidth=BAR_LW, width=0.75)
    ymax = max(counts) if counts else 1
    for rect in bars:
        h = float(rect.get_height())
        if h > 0:
            ax_top.text(rect.get_x() + rect.get_width()/2, h + max(1, 0.01*ymax),
                        f"{int(h)}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    ax_top.set_xlim(-0.5, n_cols - 0.5)
    ax_top.set_ylabel("ST counts", fontweight="bold")
    ax_top.set_xticks(range(n_cols))
    ax_top.set_xticklabels([""] * n_cols)
    ax_top.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    ax_top.yaxis.set_major_locator(MultipleLocator(500))
    for spine in ("top", "right"):
        ax_top.spines[spine].set_visible(False)
    ax_top.margins(x=0.01)

    # ----- Bottom-left: mini horizontal bars -----
    left_order = ROW_ORDER[::-1]                 # reversed then invert_yaxis() => matches matrix order
    left_y     = list(range(len(left_order)))
    totals     = [totals_by_ds.get(ds, 0) for ds in left_order]
    xmax       = max(totals) if totals else 1

    bars_left = ax_left.barh(left_y, totals, color=BAR_COLOR, edgecolor=BAR_EDGE,
                             height=0.55, linewidth=BAR_LW)

    # Bars extend left (reverse x), and flip Y to match matrix top->bottom
    ax_left.set_xlim(xmax * 1.06, 0)     # ↓ less headroom on the right
    ax_left.invert_yaxis()

    # Labels per bar (bold), smaller offset for compactness
    label_offset = 0.06 * xmax           # ↓ smaller offset
    inverted = ax_left.xaxis_inverted()
    for rect, v in zip(bars_left, totals):
        x_end    = rect.get_x() + rect.get_width()
        y_center = rect.get_y() + rect.get_height() / 2.0
        ha = 'right' if inverted else 'left'
        ax_left.text(x_end + label_offset, y_center, f"{v}", ha=ha, va="center",
                     fontsize=9, fontweight="bold", clip_on=False)

    # Cosmetic axes for left plot
    ax_left.set_yticks(left_y)
    ax_left.set_yticklabels([])
    ax_left.tick_params(axis="y", length=0, pad=4)  # ↓ slightly tighter padding
    ax_left.set_xlabel("")
    ax_left.grid(False)
    ax_left.yaxis.set_ticks_position('right')
    ax_left.yaxis.set_label_position('right')
    ax_left.spines['left'].set_visible(False)
    ax_left.spines['top'].set_visible(False)
    ax_left.spines['right'].set_visible(True)
    ax_left.margins(y=0.05)

    # ----- Bottom-right: membership matrix -----
    ax_mat.set_xlim(-0.5, n_cols - 0.5)
    ax_mat.set_ylim(-0.5, len(ROW_ORDER) - 0.5)
    ax_mat.set_yticks(range(len(ROW_ORDER)))
    ax_mat.set_yticklabels([DATASET_DISPLAY[d] for d in ROW_ORDER], fontsize=10, fontweight="bold")
    ax_mat.tick_params(axis="y", length=0, pad=8)  # ↓ tighter padding
    ax_mat.set_xticks(range(n_cols))
    ax_mat.set_xticklabels([""] * n_cols)
    ax_mat.margins(x=0.01)

    # alternating row shading
    for yi in range(len(ROW_ORDER)):
        ax_mat.axhspan(yi - 0.5, yi + 0.5, color=ROW_SHADE, alpha=0.6 if yi % 2 == 0 else 0.0, zorder=0)

    # dots & vertical connectors
    for x, sname in enumerate(x_keys):
        members = SUBSET_MEMBERS[sname]
        y_on = [i for i, ds in enumerate(ROW_ORDER) if ds in members]
        for yi in range(len(ROW_ORDER)):
            on = yi in y_on
            ax_mat.scatter(x, yi, s=DOT_SIZE,
                           facecolor=(DOT_ON_FILL if on else DOT_OFF_FILL),
                           edgecolor=DOT_EDGE, linewidth=0.6, zorder=3)
        if len(y_on) >= 2:
            ax_mat.plot([x, x], [min(y_on), max(y_on)], color=LINE_COLOR, linewidth=LINE_WIDTH, zorder=2)

    # Bold ticks everywhere
    for ax in (ax_top, ax_left, ax_mat):
        for ticklab in ax.get_xticklabels() + ax.get_yticklabels():
            ticklab.set_fontweight("bold")

    # Tighter outer margins without clipping labels
    fig.subplots_adjust(left=0.12, right=0.995, top=0.98, bottom=0.12)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ----------------------- MAIN -----------------------
def main():
    df = load_metadata(IN_PATH)
    S_A, S_B, S_C = build_sets(df)
    subsets = compute_exclusive_subsets(S_A, S_B, S_C)
    if not subsets:
        print("[WARN] All exclusive subsets are empty; nothing to plot.")
        return

    totals_by_ds = {A_LABEL: len(S_A), B_LABEL: len(S_B), C_LABEL: len(S_C)}
    plot_upset(subsets, totals_by_ds)

    print(f"[OK] Wrote: {OUT_PNG}")
    print(f"[INFO] Non-empty subsets: {list(subsets.keys())}")
    print(f"[INFO] Totals |A|={len(S_A)}, |B|={len(S_B)}, |C|={len(S_C)}, |A∪B∪C|={len(S_A|S_B|S_C)}")

if __name__ == "__main__":
    main()
