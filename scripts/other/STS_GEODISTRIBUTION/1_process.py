#!/usr/bin/env python3
"""
Generate metadata and per_ST_summary tables from KlebAMRnet, KlebNNSsero,
and klebpavia+kaspah using consistent cleaning and country→continent mapping.

Outputs:
  TASK1/results/tables/metadata.tsv
    columns: dataset, genomeID, ST, country, continent

  TASK1/results/tables/per_ST_summary.tsv
    columns: ST, n_countries, n_continents, datasets_present, subset
"""

from __future__ import annotations
from pathlib import Path
from typing import Iterable, Dict, Any, Optional, Set, List
import pandas as pd

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
TABLES_DIR = OUTPUT_DIR / "tables"

KLEBAMRNET_PATH       = FIGSHARE_DIR / "STS_GEODISTRIBUTION" / "KlebAMRnet.tsv"
KLEBPAVIA_KASPAH_PATH = FIGSHARE_DIR / "STS_GEODISTRIBUTION" / "klebpavia+kaspah.tsv"
KLEBNNSSERO_PATH      = FIGSHARE_DIR / "STS_GEODISTRIBUTION" / "KlebNNSsero.tsv"


# Dataset canonical labels (must match 'dataset' column values)
A_LABEL = "KlebAMRnet"
B_LABEL = "KlebNNSsero"
C_LABEL = "klebpavia+kaspah"

# ----------------------- NA & NORMALIZATION -----------------------
NA_TOKENS = {
    None, "", "na", "n/a", "null", "none",
    "unk", "unknown", "nan", ".", "-", "--",
    "missing", "other (international space station)"
}
NA_TOKENS = {x.lower() if isinstance(x, str) else x for x in NA_TOKENS}

def _is_na_like(x: Any) -> bool:
    if pd.isna(x):
        return True
    if isinstance(x, str):
        return x.strip().lower() in NA_TOKENS
    return False

def _norm(x: Any) -> Optional[str]:
    if _is_na_like(x):
        return None
    s = str(x).strip()
    if not s:
        return None
    return " ".join(s.split())

def drop_na_like(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    m = pd.Series(False, index=df.index)
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"Required column '{c}' missing.")
        m = m | df[c].map(_is_na_like)
    return df.loc[~m].copy()

def read_table_auto(p: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(p, sep=None, engine="python", dtype=str, keep_default_na=False)
    except Exception:
        try:
            return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        except Exception:
            return pd.read_csv(p, sep=",", dtype=str, keep_default_na=False)

# ---------------- Country → Continent mapping (full dict) ----------------
COUNTRY_TO_CONTINENT: Dict[str, str] = {
    # -------------------- OCEANIA --------------------
    "Australia": "Oceania",
    "Guam": "Oceania",
    "New Zealand": "Oceania",

    # -------------------- EUROPE --------------------
    "Albania": "Europe",
    "Austria": "Europe",
    "Belarus": "Europe",
    "Belgium": "Europe",
    "Bulgaria": "Europe",
    "Croatia": "Europe",
    "Czech Republic": "Europe",
    "Czechia": "Europe",
    "Cyprus": "Europe",
    "Denmark": "Europe",
    "England": "Europe",
    "Estonia": "Europe",
    "Finland": "Europe",
    "France": "Europe",
    "Germany": "Europe",
    "Greece": "Europe",
    "Hungary": "Europe",
    "Ireland": "Europe",
    "Italy": "Europe",
    "Latvia": "Europe",
    "Lithuania": "Europe",
    "Luxembourg": "Europe",
    "Macedonia": "Europe",  # North Macedonia
    "Malta": "Europe",
    "Montenegro": "Europe",
    "Netherlands": "Europe",
    "Norway": "Europe",
    "Poland": "Europe",
    "Portugal": "Europe",
    "Romania": "Europe",
    "Russia": "Europe",
    "Scotland": "Europe",
    "Serbia": "Europe",
    "Slovakia": "Europe",
    "Slovenia": "Europe",
    "Spain": "Europe",
    "Sweden": "Europe",
    "Switzerland": "Europe",
    "Ukraine": "Europe",
    "UK": "Europe",
    "United Kingdom": "Europe",
    "United Kingdom (England/Wales/N. Ireland)": "Europe",
    "United Kingdom (Scotland)": "Europe",
    "Wales": "Europe",

    # -------------------- ASIA --------------------
    "Afghanistan": "Asia",
    "Armenia": "Asia",
    "Bahrain": "Asia",
    "Bangladesh": "Asia",
    "Cambodia": "Asia",
    "China": "Asia",
    "Hong Kong": "Asia",
    "India": "Asia",
    "Indonesia": "Asia",
    "Iran": "Asia",
    "Iraq": "Asia",
    "Israel": "Asia",
    "Japan": "Asia",
    "Jordan": "Asia",
    "Kazakhstan": "Asia",
    "KSA": "Asia",  # alias for Saudi Arabia
    "Kuwait": "Asia",
    "Laos": "Asia",
    "Lebanon": "Asia",
    "Malaysia": "Asia",
    "Myanmar": "Asia",
    "Nepal": "Asia",
    "Oman": "Asia",
    "Pakistan": "Asia",
    "Philippines": "Asia",
    "Qatar": "Asia",
    "Saudi Arabia": "Asia",
    "Singapore": "Asia",
    "South Korea": "Asia",
    "Syria": "Asia",
    "Taiwan": "Asia",
    "Thailand": "Asia",
    "Turkey": "Asia",
    "UAE": "Asia",  # alias for United Arab Emirates
    "United Arab Emirates": "Asia",
    "Viet Nam": "Asia",
    "Vietnam": "Asia",
    "West Bank": "Asia",

    # -------------------- AFRICA --------------------
    "Algeria": "Africa",
    "Benin": "Africa",
    "Botswana": "Africa",
    "Burkina Faso": "Africa",
    "Cameroon": "Africa",
    "Djibouti": "Africa",
    "Egypt": "Africa",
    "Ethiopia": "Africa",
    "Gambia": "Africa",
    "Ghana": "Africa",
    "Kenya": "Africa",
    "Madagascar": "Africa",
    "Malawi": "Africa",
    "Morocco": "Africa",
    "Mozambique": "Africa",
    "Nigeria": "Africa",
    "Rwanda": "Africa",
    "Senegal": "Africa",
    "South Africa": "Africa",
    "Sudan": "Africa",
    "Tanzania": "Africa",
    "The Gambia": "Africa",
    "Tunisia": "Africa",
    "Uganda": "Africa",
    "Zambia": "Africa",

    # -------------------- NORTH AMERICA --------------------
    "Canada": "North America",
    "Curacao": "North America",
    "Dominican Republic": "North America",
    "Guadeloupe": "North America",
    "Guatemala": "North America",
    "Honduras": "North America",
    "Martinique": "North America",
    "Mexico": "North America",
    "United States": "North America",
    "USA": "North America",

    # -------------------- SOUTH AMERICA --------------------
    "Argentina": "South America",
    "Bolivia": "South America",
    "Brazil": "South America",
    "Chile": "South America",
    "Colombia": "South America",
    "Ecuador": "South America",
    "Paraguay": "South America",
    "Peru": "South America",
    "Uruguay": "South America",
    "Venezuela": "South America",

    # -------------------- OTHER / SPECIAL --------------------
    "Other (International Space Station)": "Other",
}

def country_to_continent(name: str) -> str:
    if not isinstance(name, str) or name.strip() == "":
        return "Other"
    return COUNTRY_TO_CONTINENT.get(name, "Other")

# ----------------------- DATASET PROCESSORS -----------------------
def process_klebamrnet(path: Path) -> pd.DataFrame:
    df = read_table_auto(path)
    rename = {"Name": "genomeID", "Country": "country", "Genotype": "ST"}
    mapped = {c: rename[c] for c in df.columns if c in rename}
    df = df.rename(columns=mapped)

    keep = [c for c in ["genomeID", "ST", "country"] if c in df.columns]
    missing = [c for c in ["genomeID", "ST", "country"] if c not in keep]
    if missing:
        raise ValueError(f"[KlebAMRnet] Missing columns after rename: {missing}")
    df = df[keep].copy()

    for c in keep:
        df[c] = df[c].map(_norm)
    df = drop_na_like(df, ["genomeID", "ST", "country"])

    df["dataset"] = A_LABEL
    return df[["dataset", "genomeID", "ST", "country"]].reset_index(drop=True)

def process_klebpavia_kaspah(path: Path) -> pd.DataFrame:
    df = read_table_auto(path)
    need = ["genomeID", "ST", "collection", "species_abbreviation"]
    keep = [c for c in need if c in df.columns]
    missing = [c for c in need if c not in keep]
    if missing:
        raise ValueError(f"[klebpavia+kaspah] Missing columns: {missing}")
    df = df[keep].copy()

    for c in keep:
        df[c] = df[c].map(_norm)

    # Restrict to K. pneumoniae (KPN)
    df = df[df["species_abbreviation"] == "KPN"].copy()

    collection_to_country = {
        "kaspah": "Australia",
        "kaspah_complete": "Australia",
        "klebpavia_published": "Italy",
        "unpublished": "Italy",
    }

    def _map_country(v: Any) -> Optional[str]:
        v = _norm(v)
        if v is None:
            return None
        key = v.replace(" ", "_").lower()
        return collection_to_country.get(key)

    df["country"] = df["collection"].map(_map_country)
    df = drop_na_like(df, ["genomeID", "ST", "country"])

    df = df[["genomeID", "ST", "country"]].copy()
    df["dataset"] = C_LABEL
    return df[["dataset", "genomeID", "ST", "country"]].reset_index(drop=True)

def process_klebnnssero(path: Path) -> pd.DataFrame:
    df = read_table_auto(path)
    rename = {"Genome Name": "genomeID", "Country": "country", "ST": "ST"}
    mapped = {c: rename[c] for c in df.columns if c in rename}
    df = df.rename(columns=mapped)

    keep = [c for c in ["genomeID", "ST", "country"] if c in df.columns]
    missing = [c for c in ["genomeID", "ST", "country"] if c not in keep]
    if missing:
        raise ValueError(f"[KlebNNSsero] Missing columns after rename: {missing}")
    df = df[keep].copy()

    for c in keep:
        df[c] = df[c].map(_norm)
    df = drop_na_like(df, ["genomeID", "ST", "country"])

    df["dataset"] = B_LABEL
    return df[["dataset", "genomeID", "ST", "country"]].reset_index(drop=True)

# ----------------------- MERGE + ENRICH -----------------------
def load_clean_merge_add_continent() -> pd.DataFrame:
    print("[1/3] Processing KlebAMRnet…")
    a = process_klebamrnet(KLEBAMRNET_PATH)

    print("[2/3] Processing klebpavia+kaspah…")
    c = process_klebpavia_kaspah(KLEBPAVIA_KASPAH_PATH)

    print("[3/3] Processing KlebNNSsero…")
    b = process_klebnnssero(KLEBNNSSERO_PATH)

    combined = pd.concat([a, b, c], ignore_index=True)

    # Canonicalize strings & drop duplicates (final safety)
    for col in ["dataset", "genomeID", "ST", "country"]:
        combined[col] = combined[col].astype(str).map(lambda s: " ".join(s.split()))
    combined = combined.drop_duplicates(subset=["dataset", "genomeID", "ST", "country"]).reset_index(drop=True)

    # Add continent
    combined["continent"] = combined["country"].map(country_to_continent)

    return combined

# ----------------------- PER-ST SUMMARY -----------------------
def compute_per_st_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns one row per ST with:
      - n_countries, n_continents
      - datasets_present (pipe-separated dataset labels)
      - subset (A unique, B unique, C unique, A∩B only, A∩C only, B∩C only, A∩B∩C)
    """
    n_countries = df.groupby("ST")["country"].nunique()
    n_continents = df.groupby("ST")["continent"].nunique()

    present_by_st = df.groupby("ST")["dataset"].apply(lambda s: set(s.unique()))

    def subset_label(present: Set[str]) -> str:
        hasA = A_LABEL in present
        hasB = B_LABEL in present
        hasC = C_LABEL in present
        if hasA and hasB and hasC:
            return "A∩B∩C"
        if hasA and hasB and not hasC:
            return "A∩B only"
        if hasA and hasC and not hasB:
            return "A∩C only"
        if hasB and hasC and not hasA:
            return "B∩C only"
        if hasA and not hasB and not hasC:
            return "A unique"
        if hasB and not hasA and not hasC:
            return "B unique"
        if hasC and not hasA and not hasB:
            return "C unique"
        return "Other"

    rows: List[Dict[str, Any]] = []
    for st, present in present_by_st.items():
        rows.append({
            "ST": st,
            "n_countries": int(n_countries.get(st, 0)),
            "n_continents": int(n_continents.get(st, 0)),
            "datasets_present": "|".join(sorted(present)),
            "subset": subset_label(present),
        })

    per_st = pd.DataFrame(rows, columns=["ST", "n_countries", "n_continents", "datasets_present", "subset"])
    return per_st.sort_values(["subset", "ST"]).reset_index(drop=True)

# ----------------------- MAIN -----------------------
def main():
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    # 1) Load, clean, merge, add continent
    combined = load_clean_merge_add_continent()

    # 2) metadata.tsv
    meta_out = TABLES_DIR / "metadata.tsv"
    combined[["dataset", "genomeID", "ST", "country", "continent"]].to_csv(meta_out, sep="\t", index=False)
    print(f"[OK] Wrote metadata: {meta_out} ({len(combined):,} rows)")

    # 3) per_ST_summary.tsv
    per_st = compute_per_st_summary(combined)
    per_st_out = TABLES_DIR / "per_ST_summary.tsv"
    per_st.to_csv(per_st_out, sep="\t", index=False)
    print(f"[OK] Wrote per_ST_summary: {per_st_out} ({len(per_st):,} STs)")

if __name__ == "__main__":
    main()
