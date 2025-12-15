#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------- YAML LOADER -----------------------
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

# ----------------------- PATHS -----------------------
IN_PATH  = Path(TABLES_DIR, "metadata.tsv")  # must have: dataset, genomeID, ST, country, continent
OUT_PNG  = Path(PLOTS_DIR, "stacked_barplot.png")

# ----------------------- DATASET & LABELS -----------------------
A_LABEL = "KlebAMRnet"        # Pathogenwatch
B_LABEL = "KlebNNSsero"
C_LABEL = "klebpavia+kaspah"  # GWAS

SUBSET_ORDER = ["A∩B∩C", "A∩C only", "A∩B only", "A unique", "B unique", "C unique"]
SUBSET_LABEL_MAP: Dict[str,str] = {
    "A∩B∩C":    "KlebAMRnet\n∩\nKlebNNSsero\n∩\nGWAS",
    "A∩C only": "KlebAMRnet\n∩\nGWAS",
    "A∩B only": "KlebAMRnet\n∩\nKlebNNSsero",
    "A unique": "KlebAMRnet",
    "B unique": "KlebNNSsero",
    "C unique": "GWAS",
}

# ----------------------- BINS & COLORS -----------------------
COUNTRY_BINS   = ["1", "2–3", "4–9", "≥10"]
COUNTRY_DRAW   = ["≥10", "4–9", "2–3", "1"]  # draw from bottom upwards
COUNTRY_COLORS = {"≥10":"#B22222","4–9":"#E2583E","2–3":"#F7CDBB","1":"#D9D9D9"}

CONTINENT_BINS   = ["1", "2–3", "4–5", "≥6"]
CONTINENT_DRAW   = ["≥6", "4–5", "2–3", "1"]
CONTINENT_COLORS = {"≥6":"#0B4DA2","4–5":"#2F7DC1","2–3":"#79B1D9","1":"#D9D9D9"}

FIGSIZE = (8, 6)
Y_TICKS = [0.0, 0.25, 0.5, 0.75, 1.0]

# ----------------------- IO PREP -----------------------
PLOTS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------- HELPERS -----------------------
def load_metadata(p: Path) -> pd.DataFrame:
    df = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
    need = {"dataset","genomeID","ST","country","continent"}
    miss = need - set(df.columns)
    if miss:
        raise ValueError(f"Missing required columns: {sorted(miss)}")
    return df

def build_sets(df: pd.DataFrame) -> Tuple[Set[str],Set[str],Set[str]]:
    return (
        set(df.loc[df["dataset"]==A_LABEL, "ST"].unique()),
        set(df.loc[df["dataset"]==B_LABEL, "ST"].unique()),
        set(df.loc[df["dataset"]==C_LABEL, "ST"].unique()),
    )

def breadth_by(df: pd.DataFrame, st_subset: Set[str], scope: Set[str], field: str) -> pd.Series:
    if not st_subset: return pd.Series(dtype=int)
    scoped = df[df["dataset"].isin(scope) & df["ST"].isin(st_subset)][["ST", field]]
    if scoped.empty: return pd.Series(dtype=int)
    return scoped.groupby("ST")[field].nunique().astype(int)

def bin_countries(n: int) -> str:
    return "1" if n<=1 else "2–3" if n<=3 else "4–9" if n<=9 else "≥10"

def bin_continents(n: int) -> str:
    return "1" if n<=1 else "2–3" if n<=3 else "4–5" if n<=5 else "≥6"

# ----------------------- MAIN -----------------------
def main():
    df = load_metadata(IN_PATH)

    S_A, S_B, S_C = build_sets(df)
    subset_defs = [
        ("A∩B∩C",   (S_A & S_B & S_C),          {A_LABEL, B_LABEL, C_LABEL}),
        ("A∩C only",((S_A & S_C) - S_B),        {A_LABEL, C_LABEL}),
        ("A∩B only",((S_A & S_B) - S_C),        {A_LABEL, B_LABEL}),
        ("A unique",(S_A - (S_B | S_C)),        {A_LABEL}),
        ("B unique",(S_B - (S_A | S_C)),        {B_LABEL}),
        ("C unique",(S_C - (S_A | S_B)),        {C_LABEL}),
    ]

    recs = []
    for label, st_set, scope in subset_defs:
        if not st_set: continue
        # Countries
        bc = breadth_by(df, st_set, scope, "country")
        if not bc.empty:
            bins = bc.map(bin_countries); counts = bins.value_counts().to_dict(); total = int(bc.shape[0])
            for b in COUNTRY_BINS:
                recs.append({"subset":label,"level":"countries","bin":b,"count_STs":int(counts.get(b,0)),"n_total_STs":total})
        # Continents
        bt = breadth_by(df, st_set, scope, "continent")
        if not bt.empty:
            bins = bt.map(bin_continents); counts = bt.value_counts().to_dict(); total = int(bt.shape[0])  # typo fixed below
            # ^ OOPS — fix line:
            bins = bt.map(bin_continents); counts = bins.value_counts().to_dict(); total = int(bt.shape[0])
            for b in CONTINENT_BINS:
                recs.append({"subset":label,"level":"continents","bin":b,"count_STs":int(counts.get(b,0)),"n_total_STs":total})

    if not recs:
        print("[WARN] No data to plot.")
        return

    out = pd.DataFrame.from_records(recs)
    out["proportion"] = out["count_STs"] / out["n_total_STs"].where(out["n_total_STs"]>0, 1)

    present = [s for s in SUBSET_ORDER if s in set(out["subset"])]
    x_labels = [SUBSET_LABEL_MAP.get(s, s) for s in present]

    def pivot_for(level: str, bins: List[str]):
        df_lv = out[(out["level"]==level) & (out["subset"].isin(present)) & (out["bin"].isin(bins))]
        idx = pd.MultiIndex.from_product([present, bins], names=["subset","bin"])
        df_lv = df_lv.set_index(["subset","bin"]).reindex(idx, fill_value=0.0)
        tbl = df_lv["proportion"].unstack("bin").loc[present, bins]
        totals = (df_lv.reset_index().drop_duplicates(["subset"])
                  .set_index("subset")["n_total_STs"].reindex(present).astype(int).tolist())
        return tbl, totals

    P_countries, totals_c = pivot_for("countries", COUNTRY_BINS)
    P_continents, _      = pivot_for("continents", CONTINENT_BINS)

    # ---- Figure (no legend column), zero vertical gap ----
    fig = plt.figure(figsize=FIGSIZE)
    gs  = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1,1], hspace=0.1)

    ax_top = fig.add_subplot(gs[0,0])
    ax_bot = fig.add_subplot(gs[1,0])

    # ---------------- TOP (countries) ----------------
    bottoms = [0.0]*len(present)
    for b in COUNTRY_DRAW:
        vals = P_countries[b].tolist()
        ax_top.bar(x_labels, vals, bottom=bottoms,
                   color=COUNTRY_COLORS[b], edgecolor="black", linewidth=0.6, label=b)
        bottoms = [a+v for a,v in zip(bottoms, vals)]
    ax_top.set_ylim(0,1.06); ax_top.set_yticks(Y_TICKS)
    ax_top.grid(axis="y", color='black', linestyle=":", linewidth=0.6, alpha=0.6)

    # Top: hide x-ticks/labels
    ax_top.set_xticks([])
    ax_top.tick_params(axis="x", which="both", length=0, labelbottom=False)

    # Totals on top of bars (smaller font)
    for x, n in enumerate(totals_c):
        ax_top.text(x, 1.1, f"n={n}", ha="center", va="bottom",
                    fontsize=10, fontweight="bold")  # smaller label here

    # ---------------- BOTTOM (continents) ----------------
    bottoms = [0.0]*len(present)
    for b in CONTINENT_DRAW:
        vals = P_continents[b].tolist()
        ax_bot.bar(x_labels, vals, bottom=bottoms,
                   color=CONTINENT_COLORS[b], edgecolor="black", linewidth=0.6, label=b)
        bottoms = [a+v for a,v in zip(bottoms, vals)]
    ax_bot.set_ylim(0,1.06); ax_bot.set_yticks(Y_TICKS)
    ax_bot.grid(axis="y", color='black', linestyle=":", linewidth=0.6, alpha=0.6)

    # ---------------- Shared formatting ----------------
    for ax in (ax_top, ax_bot):
        ax.set_xlabel("")
        ax.set_title("")
        ax.tick_params(axis="x", labelsize=9)
        # Bold all tick labels on both axes
        for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
            lab.set_fontweight("bold")

    # Shared Y label
    fig.text(0, 0.5, "Proportion of STs", va="center", rotation="vertical", fontweight="bold")

    # Tight, with zero vertical gap
    fig.subplots_adjust(left=0.10, right=0.98, top=0.98, bottom=0.12, hspace=1)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print(f"[OK] Wrote: {OUT_PNG}")

if __name__ == "__main__":
    main()
