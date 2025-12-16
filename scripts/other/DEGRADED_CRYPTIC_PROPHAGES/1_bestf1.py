# Combine HQ and LQ pyseer hits into a single long-form table for lasso-mode K-locus F1 comparison

import pandas as pd
from pathlib import Path
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
CONFIG_PATH = Path('../../../config/config.yaml')
CONFIG_DICT = load_yaml_dict(CONFIG_PATH)

FIGSHARE_DATA_DIR = Path(CONFIG_DICT['paths']['janusz']['main']) / Path(CONFIG_DICT['paths']['janusz']['figshare_data'])
FIGSHARE_GWAS_DIR = Path(CONFIG_DICT['paths']['janusz']['main']) / Path(CONFIG_DICT['paths']['janusz']['figshare_gwas'])
OUTPUT_DIR = Path(CONFIG_DICT['paths']['janusz']['main']) / Path(CONFIG_DICT['paths']['janusz']['output'])

OUTPUT_DIR = OUTPUT_DIR / "DEGRADED_CRYPTIC_PROPHAGES"
FIGURES_DIR = OUTPUT_DIR / "plots"
TABLES_DIR = OUTPUT_DIR / "tables"

# Input files 
hq_path = Path(FIGSHARE_DATA_DIR, 'REVIEW/DEGRADED_CRYPTIC_PROPHAGES/GWAS_HIGH_QUALITY_PROPHAGES/3_PROCESSING/pyseer_hits_all.tsv')
lq_path = Path(FIGSHARE_GWAS_DIR, 'GWAS/3_PROCESSING/pyseer_hits_all.tsv')

out_path = Path(TABLES_DIR, 'bestf1.tsv')

Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)
Path(TABLES_DIR).mkdir(parents=True, exist_ok=True)

# load
hq = pd.read_csv(hq_path, sep="\t", dtype=str)
lq = pd.read_csv(lq_path, sep="\t", dtype=str)

hq["quality"] = "hq"; lq["quality"] = "lq"
df_all = pd.concat([hq, lq], ignore_index=True)
df = df_all[df_all.get("mode", "") == "lasso"].copy()

# Ensure helper numeric columns
def to_num(s):
    return pd.to_numeric(s, errors="coerce")

if "F1_score" in df.columns and "F1_score_num" not in df.columns:
    df["F1_score_num"] = to_num(df["F1_score"])

# Optional secondary sort key if present
if "minus_log10_pvalue_corr" in df.columns:
    df["_minuslog_num"] = to_num(df["minus_log10_pvalue_corr"])
else:
    df["_minuslog_num"] = pd.NA

# Keep only rows with non-null F1 for the max selection
df_nonan = df[df["F1_score_num"].notna()].copy()

# Sorting to impose deterministic tie-breaks before groupby.head(1)
sort_keys = ["version", "locus", "quality", "F1_score_num", "_minuslog_num"]
ascending = [True, True, True, False, False]
present_keys = [k for k in sort_keys if k in df_nonan.columns]
asc = [ascending[sort_keys.index(k)] for k in present_keys]

df_sorted = df_nonan.sort_values(by=present_keys, ascending=asc, kind="mergesort")

group_keys = [k for k in ["version", "locus", "quality"] if k in df_sorted.columns]

best_rows = (
    df_sorted.groupby(group_keys, as_index=False, sort=False)
    .head(1)
    .reset_index(drop=True)
)

# Clean up helper (keep all original columns otherwise)
if "_minuslog_num" in best_rows.columns:
    best_rows = best_rows.drop(columns=["_minuslog_num"])

# ---------- COMPLETION STEP (simple) ----------
# Ensure every (locus, quality) has all six clustering levels; add empty rows where missing
required_versions = ["PCI00C50", "PCI00C80", "PCI50C50", "PCI50C80", "PCI80C50", "PCI80C80"]

# All observed (locus, quality) pairs in the lasso source
pairs = df[["locus", "quality"]].dropna().drop_duplicates()

# Build full grid via a simple cartesian merge (no fancy functions)
pairs["_k"] = 1
vers_df = pd.DataFrame({"version": required_versions})
vers_df["_k"] = 1
grid = pairs.merge(vers_df, on="_k").drop(columns="_k")[["version", "locus", "quality"]]

# Left-join grid with best_rows to add empty rows where no PC was selected
best_rows = grid.merge(best_rows, on=["version", "locus", "quality"], how="left")

# ------------------------------------------------

# Save
best_rows.to_csv(out_path, sep="\t", index=False)

# Provide a short summary as structured output
{
    "groups": len(best_rows),
    "columns": list(best_rows.columns),
    "output_file": str(out_path)
}
