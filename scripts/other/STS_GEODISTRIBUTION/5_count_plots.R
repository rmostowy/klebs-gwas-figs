#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

# ---- Locate config.yaml (../../../config/config.yaml) ----
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  getwd()
}
SCRIPT_DIR  <- get_script_dir()
CONFIG_PATH <- normalizePath(file.path(SCRIPT_DIR, "..", "..", "..", "config", "config.yaml"), mustWork = TRUE)

# ---- Read config and construct paths ----
cfg          <- yaml::read_yaml(CONFIG_PATH)
FIGSHARE_DIR <- file.path(cfg$paths$janusz$main, cfg$paths$janusz$figshare_data, "REVIEW")
OUTPUT_DIR   <- file.path(cfg$paths$janusz$main, cfg$paths$janusz$output, "STS_GEODISTRIBUTION")
TABLES_DIR   <- file.path(OUTPUT_DIR, "tables")
PLOTS_DIR    <- file.path(OUTPUT_DIR, "plots")
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Inputs ----
metadata.filename <- file.path(TABLES_DIR, "metadata.tsv")
if (!file.exists(metadata.filename)) stop(sprintf("metadata.tsv not found at: %s", metadata.filename))
message("[INFO] Using metadata: ", metadata.filename)
message("[INFO] Plots will be written to: ", PLOTS_DIR)

# ---- Load and check columns ----
metadata <- fread(metadata.filename)
needed <- c("dataset", "ST", "country", "continent")
miss <- setdiff(needed, names(metadata))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

# ---- Dataset labels ----
DS_GWAS <- "klebpavia+kaspah"

# ---- Compute bins (excluding GWAS countries/continents where appropriate) ----
sts.gwas        <- unique(metadata[dataset == DS_GWAS, ST])
countries.gwas  <- unique(metadata[dataset == DS_GWAS, country])
continents.gwas <- unique(metadata[dataset == DS_GWAS, continent])

gwas.st.distr.list <- lapply(seq_along(sts.gwas), function(st.index){
  this.st <- sts.gwas[st.index]
  this.st.metadata <- metadata[ST %in% this.st]

  this.st.metadata.gwas.countr.excl <- this.st.metadata[!country %in% countries.gwas]
  this.st.metadata.gwas.cont.excl   <- this.st.metadata[!continent %in% continents.gwas]

  no.countries                     <- uniqueN(this.st.metadata$country)
  no.continents                    <- uniqueN(this.st.metadata$continent)
  no.countries.gwas.countr.excl    <- uniqueN(this.st.metadata.gwas.countr.excl$country)
  no.continents.gwas.countr.excl   <- uniqueN(this.st.metadata.gwas.countr.excl$continent)
  no.countries.gwas.cont.excl      <- uniqueN(this.st.metadata.gwas.cont.excl$country)
  no.continents.gwas.cont.excl     <- uniqueN(this.st.metadata.gwas.cont.excl$continent)

  data.table(no.countries, no.continents,
             no.countries.gwas.countr.excl, no.continents.gwas.countr.excl,
             no.countries.gwas.cont.excl,   no.continents.gwas.cont.excl)
})
gwas.st.distr <- rbindlist(gwas.st.distr.list)

gwas.st.distr[, n_countries_bin := fifelse(
  no.countries.gwas.countr.excl >= 6, "6+", as.character(no.countries.gwas.countr.excl)
)]
gwas.st.distr[, n_continents_bin := fifelse(
  no.continents.gwas.cont.excl >= 2, "2+", as.character(no.continents.gwas.cont.excl)
)]

gwas.st.distr[, n_countries_bin := factor(n_countries_bin, levels = c("0","1","2","3","4","5","6+"))]
gwas.st.distr[, n_continents_bin := factor(n_continents_bin, levels = c("0","1","2+"))]

dt_plot_countries  <- gwas.st.distr[, .N, by = n_countries_bin]
dt_plot_continents <- gwas.st.distr[, .N, by = n_continents_bin]

# ---- Bold theme + anti-clipping tweaks ----
bold_theme <- theme_classic(base_size = 10) +
  theme(
    text              = element_text(face = "bold"),
    axis.title        = element_text(face = "bold"),
    axis.text         = element_text(face = "bold"),
    plot.title        = element_text(face = "bold", hjust = 0.5, margin = margin(b = 8)),
    legend.title      = element_text(face = "bold"),
    legend.text       = element_text(face = "bold"),
    strip.text        = element_text(face = "bold"),
    plot.margin       = margin(t = 10, r = 18, b = 10, l = 12)  # more room so labels don't get cut
  )

# Extra headroom so tall labels don’t touch the edge; clip off so annotations can hang outside
y_expand <- expansion(mult = c(0, 0.08))

p.countries <- ggplot(dt_plot_countries, aes(x = n_countries_bin, y = N)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(expand = y_expand) +
  labs(
    x = "Number of countries (excluding Italy and Australia)",
    y = "Number of sequence types (STs)",
    title = "Geographic spread of STs outside Italy and Australia"
  ) +
  coord_cartesian(clip = "off") +
  bold_theme

p.continents <- ggplot(dt_plot_continents, aes(x = n_continents_bin, y = N)) +
  geom_col(fill = "steelblue") +
  scale_y_continuous(expand = y_expand) +
  labs(
    x = "Number of continents\n(excluding Europe and Oceania)",
    y = "Number of sequence types (STs)",
    title = "Geographic spread of STs\noutside Europe and Oceania"
  ) +
  coord_cartesian(clip = "off") +
  bold_theme

# ---- Save (cairo_pdf to embed fonts nicely) ----
ggsave(file.path(PLOTS_DIR, "plot1_countries.pdf"),  p.countries,  width = 6, height = 4)
ggsave(file.path(PLOTS_DIR, "plot2_continents.pdf"), p.continents, width = 3, height = 4)

message("[OK] Wrote: ", file.path(PLOTS_DIR, "plot1_countries.pdf"))
message("[OK] Wrote: ", file.path(PLOTS_DIR, "plot2_continents.pdf"))
