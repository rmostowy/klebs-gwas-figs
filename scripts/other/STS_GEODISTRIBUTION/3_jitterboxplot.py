#!/usr/bin/env python3
"""
Jitter + box plot (pooled GWAS subsets):
X = exact continents per ST (1,2,3,4,5,6; with 6 == >=6)
Y = countries per ST

Filters STs to the union of:
  - KlebAMRnet|KlebNNSsero|klebpavia+kaspah
  - KlebAMRnet|klebpavia+kaspah
  - klebpavia+kaspah

Reads:
  /Users/januszkoszucki/MGG Dropbox/Projects/kleb-prophage-div/2025-02-12_KLEBDATA_LIGHT/3_GWAS/REVIEW/TASK1/results/tables/per_ST_summary.tsv

Writes:
  /Users/januszkoszucki/MGG Dropbox/Projects/kleb-prophage-div/2025-02-12_KLEBDATA_LIGHT/3_GWAS/REVIEW/TASK1/results/plots/jitter_box_countries_by_continents_selected_gwas.pdf
"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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
CONFIG_PATH = Path('../../../config/config.yaml')
CONFIG_DICT = load_yaml_dict(CONFIG_PATH)

FIGSHARE_DIR = Path(CONFIG_DICT['paths']['janusz']['main']) / Path(CONFIG_DICT['paths']['janusz']['figshare_data']) / 'REVIEW'
OUTPUT_DIR = Path(CONFIG_DICT['paths']['janusz']['main']) / Path(CONFIG_DICT['paths']['janusz']['output'])
OUTPUT_DIR = OUTPUT_DIR / "STS_GEODISTRIBUTION"
PLOTS_DIR = OUTPUT_DIR / "plots"
TABLES_DIR = OUTPUT_DIR / "tables"


# ----------------------- PATHS -----------------------
INPUT_TABLE = TABLES_DIR / "per_ST_summary.tsv"
OUT_FILE    = PLOTS_DIR / "jitter_box.pdf"


# ----------------------- STYLE -----------------------
def set_plot_style() -> None:
    plt.rcParams.update({
        "figure.figsize": (10, 6),
        "savefig.dpi": 300,
        "font.size": 11,
        "font.weight": "bold",
        "axes.labelsize": 14,
        "axes.labelweight": "bold",
        "axes.titlesize": 14,
        "axes.titleweight": "bold",
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 10,
    })

# make dir
Path(PLOTS_DIR).mkdir(exist_ok=True, parents=True)
Path(OUTPUT_DIR).mkdir(exist_ok=True, parents=True)

def bold_ticklabels(ax: plt.Axes) -> None:
    for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        try:
            lab.set_fontweight("bold")
        except Exception:
            pass

# ----------------------- CONFIG -----------------------
SUBSETS = {
    "KlebAMRnet|KlebNNSsero|klebpavia+kaspah",
    "KlebAMRnet|klebpavia+kaspah",
    # "klebpavia+kaspah",
}

# Palette mapped by exact count: 1, 2–3, 4–5, 6(≥6)
CONTINENT_COLORS = {"≥6": "#08519C", "4–5": "#3182BD", "2–3": "#6BAED6", "1": "#D9D9D9"}
def color_for_count(n_continents_capped: int) -> str:
    if n_continents_capped == 1: return CONTINENT_COLORS["1"]
    if n_continents_capped in (2, 3): return CONTINENT_COLORS["2–3"]
    if n_continents_capped in (4, 5): return CONTINENT_COLORS["4–5"]
    return CONTINENT_COLORS["≥6"]  # 6 or more

RNG_SEED = 42
JITTER_WIDTH = 0.18
ALPHA = 0.75   # lower alpha to separate hues better

# ----------------------- IO -----------------------
def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input table not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype={"ST": str, "datasets_present": str}, keep_default_na=False)
    df["n_countries"]  = pd.to_numeric(df.get("n_countries"), errors="coerce")
    df["n_continents"] = pd.to_numeric(df.get("n_continents"), errors="coerce")

    # Select GWAS-related union; drop missing; cast to ints
    df = df[df["datasets_present"].isin(SUBSETS)].dropna(subset=["n_countries", "n_continents"]).copy()
    df["n_countries"]  = df["n_countries"].astype(int)
    df["n_continents"] = df["n_continents"].astype(int)

    # Cap continents at 6 (6 == >=6)
    df["n_continents_capped"] = df["n_continents"].clip(upper=6)
    return df

# ----------------------- PLOT -----------------------
def make_jitter_box_plot(df: pd.DataFrame, out_path: Path) -> None:
    if df.empty:
        print("[WARN] No data after filtering; nothing to plot.")
        return

    set_plot_style()
    fig, ax = plt.subplots()

    # X groups are exact counts 1..6 (only those present)
    x_levels_all = [1, 2, 3, 4, 5, 6]
    present_levels = [k for k in x_levels_all if (df["n_continents_capped"] == k).any()]
    pos = {k: i for i, k in enumerate(present_levels, start=1)}

    # --- Boxplot: unfilled black, alpha reduced ---
    data_for_box = [df.loc[df["n_continents_capped"] == k, "n_countries"].values for k in present_levels]
    bp = ax.boxplot(
        data_for_box,
        positions=[pos[k] for k in present_levels],
        widths=0.5,
        patch_artist=True,
        showmeans=False,
        medianprops={"linewidth": 2, "color": "black", "alpha": ALPHA},
        boxprops={"linewidth": 1.2, "edgecolor": "black"},
        whiskerprops={"linewidth": 1.2, "color": "black", "alpha": ALPHA},
        capprops={"linewidth": 1.2, "color": "black", "alpha": ALPHA},
        flierprops={"marker": "", "markersize": 0},
    )
    for patch in bp["boxes"]:
        patch.set_facecolor("none")
        patch.set_edgecolor("black")
        patch.set_alpha(ALPHA)

    # --- Jittered points: colored by exact count group, alpha reduced ---
    rng = np.random.default_rng(RNG_SEED)
    for k in present_levels:
        sub = df[df["n_continents_capped"] == k]
        y_vals = sub["n_countries"].to_numpy()
        if y_vals.size == 0:
            continue
        x0 = pos[k]
        x_jit = x0 + rng.uniform(-JITTER_WIDTH, JITTER_WIDTH, size=y_vals.size)
        ax.scatter(
            x_jit, y_vals,
            s=40,
            alpha=ALPHA,
            edgecolor="black",
            linewidth=0.5,
            c=color_for_count(k),
        )

    # Labels, ticks, limits
    ax.set_xlabel("\nDistinct continents per ST", fontweight="bold")
    ax.set_ylabel("\nDistinct countries per ST", fontweight="bold")
    # ax.set_title("Countries vs Continents per ST (GWAS pooled with overlaps)", fontweight="bold")

    ax.set_xticks([pos[k] for k in present_levels])

    # --- Add ST counts under each X tick: "k\n(n=NN)" ---
    counts_per_level = df["n_continents_capped"].value_counts().to_dict()
    xtick_labels = [f"{k}\n(n={counts_per_level.get(k, 0)})" for k in present_levels]
    ax.set_xticklabels(xtick_labels)

    ax.set_xlim(0.5, len(present_levels) + 0.5)
    ax.set_ylim(bottom=0)

    # Grid & aesthetics
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    bold_ticklabels(ax)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Wrote: {out_path}")

def main():
    df = load_table(INPUT_TABLE)
    make_jitter_box_plot(df, OUT_FILE)

if __name__ == "__main__":
    main()
