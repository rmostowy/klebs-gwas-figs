### load libraries
suppressWarnings(library(yaml))
suppressWarnings(library(data.table))
suppressWarnings(library(ape))
suppressWarnings(library(ggtree))
suppressWarnings(library(ggtreeExtra))
suppressWarnings(library(ggplot2))
suppressWarnings(library(RColorBrewer))

### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.yaml")
cfg <- yaml::read_yaml(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

### LOAD DATA ###

data.path <- file.path(cfg$paths$rafal$main, cfg$paths$rafal$db_output_rel)
bacteria.metadata.filename <- file.path(data.path, cfg$data$bact_metadata)
bacteria.metadata <- fread(bacteria.metadata.filename)

prophage.metadata.filename <- file.path(data.path, cfg$data$prophage_metadata)
prophage.metadata <- fread(prophage.metadata.filename)

### LOAD PARMS ###

SC.CUTOFF <- cfg$params$sc_cutoff
SC.CUTOFF.GWAS <- cfg$params$sc_cutoff_gwas

### PROCESS DATA ###

# index
setkey(bacteria.metadata, genomeID)

# Fix collection information
bacteria.metadata[collection == "kaspah"]$collection <- "KASPAH"
bacteria.metadata[collection == "kaspah_complete"]$collection <- "KASPAH-REF"
bacteria.metadata[collection == "klebpavia_published"]$collection <- "KLEBPAVIA"
bacteria.metadata[collection == "unpublished"]$collection <- "KLEBPAVIA"

# Calculate bacterial diversity per K locus
k_locus_diversity <- bacteria.metadata[, .(total_SCs = uniqueN(MGG_SC)), by = MGG_K_locus]
k_locus_kaspah_ref <- bacteria.metadata[collection == "KASPAH-REF", .(KASPAH_REF_SCs = uniqueN(MGG_SC)), by = MGG_K_locus]
k_locus_summary <- merge(k_locus_diversity, k_locus_kaspah_ref, by = "MGG_K_locus", all.x = TRUE)
k_locus_summary <- k_locus_summary[MGG_K_locus != "UNK"]
k_locus_summary[is.na(KASPAH_REF_SCs), KASPAH_REF_SCs := 0]
k_locus_summary[, other_SCs := total_SCs - KASPAH_REF_SCs]

# Calculate prophage diversity
phage_summary <- prophage.metadata[, .(
  no.pv.all = uniqueN(wgrr95)
), by = MGG_K_locus]
phage_summary_kaspah_ref <- prophage.metadata[genomeID %in% bacteria.metadata[collection == "KASPAH-REF", genomeID], .(
  no.pv.kaspah.ref = uniqueN(wgrr95)
), by = MGG_K_locus]
phage_summary <- merge(phage_summary, phage_summary_kaspah_ref, by = "MGG_K_locus", all.x = TRUE)
phage_summary[is.na(no.pv.all), no.pv.all := 0]
phage_summary[is.na(no.pv.kaspah.ref), no.pv.kaspah.ref := 0]
phage_summary[, no.pv.other := no.pv.all - no.pv.kaspah.ref]
phage_summary <- phage_summary[MGG_K_locus != "UNK"]

# Merge diversity data
diversity.data <- merge(k_locus_summary, phage_summary, by = "MGG_K_locus")
diversity.data <- diversity.data[total_SCs >= SC.CUTOFF]
k.locus.order <- diversity.data[order(-total_SCs)]$MGG_K_locus

# Convert to long format for plotting both SCs and Phage Variants correctly
diversity.data.melted <- melt(diversity.data, 
                              id.vars = "MGG_K_locus", 
                              measure.vars = c("KASPAH_REF_SCs", "other_SCs", 
                                               "no.pv.kaspah.ref", "no.pv.other"),
                              variable.name = "category",
                              value.name = "count")

# Create a new column to distinguish between SCs and Phage Variants
diversity.data.melted[, facet_group := ifelse(category %in% c("KASPAH_REF_SCs", "other_SCs"), 
                                              "Sequence Clusters", "Phage Variants")]

# Standardize the SC source labels
diversity.data.melted[, SC_source := ifelse(category %in% c("KASPAH_REF_SCs", "no.pv.kaspah.ref"), 
                                            "KASPAH-REF", "Other")]

###############################
## Plot Figure 1B
###############################

col.palette <- brewer.pal(12, "Paired")
col.blue <- col.palette[1]
col.red <- col.palette[6]
col.brown <- col.palette[12]

data.colors <- c("KASPAH-REF" = col.red,  # Dark red
                 "Other" = col.blue)       # Dark green
facet.order <- c("Sequence Clusters", "Phage Variants")

diversity.data.melted$MGG_K_locus <- factor(diversity.data.melted$MGG_K_locus, levels = k.locus.order)
diversity.data.melted$facet_group <- factor(diversity.data.melted$facet_group, levels = facet.order)
diversity.data.melted[, SC_source := factor(SC_source, levels = c("Other", "KASPAH-REF"))]

p.panelB <- ggplot(diversity.data.melted, aes(x = MGG_K_locus, y = count, fill = SC_source)) +
  geom_bar(stat = "identity", color = "black", position = "stack", width = 0.7) +
  facet_grid(facet_group ~ ., scales = "free_y") +  # Separate SCs and PVs into different facets
  scale_fill_manual(values = data.colors) +  # Red & Green
  theme_minimal(base_size = 14) +
  # ✅ Add geom_hline() but restrict it to "Sequence Clusters" facet
  geom_hline(data = subset(diversity.data.melted, facet_group == "Sequence Clusters"), 
             aes(yintercept = SC.CUTOFF.GWAS), 
             linetype = "dashed", color = col.brown, linewidth = 1.2) +
  labs(
    x = "K-locus",
    y = "Count",
    fill = "SC Source"
  ) +
  theme(
    text = element_text(family = "Myriad Pro"),
    panel.grid.major.x = element_line(color = "gray80", linetype = "dotted"),
    panel.grid.major.y = element_line(color = "gray80"), 
    strip.text = element_text(face = "bold", size = 18),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 13),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )
print(p.panelB)

# ggsave("Figure_1B.png", p.panelB, width = 11, height = 8, device = cairo_pdf)
ggsave("Figure_1B.png", p.panelB, width = 11, height = 8)

###############################
## Plot Supplementary Figure
###############################

dir.create("supplementary figures")

# Load the K-locus Summary Table (containing number of SCs per K-locus)
# Assuming the table name is `k_locus_summary`
k_locus_summary <- k_locus_summary[total_SCs >= SC.CUTOFF]

# Merge completeness data with K-locus summary to ensure proper sorting
prophage_completeness <- merge(prophages[, .(MGG_K_locus, completeness)], k_locus_summary, by = "MGG_K_locus")

# Plot the distribution of prophage completeness per K-locus (Jitter Plot)
p.completeness.all <- ggplot(prophage_completeness, aes(x = reorder(MGG_K_locus, -total_SCs), y = completeness)) +
  geom_boxplot(outlier.shape = NA, fill = "#A6CEE3", color = "black", alpha = 0.6) + 
  geom_jitter(width = 0.2, alpha = 0.5, color = "#1A85FF") +  # Blue points for clarity
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Prophage Completeness Across K-loci",
    x = "K-locus",
    y = "Prophage Completeness (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Compute phage variant completeness by selecting the max completeness for each unique phage variant
phage_variant_completeness <- prophages[, .(phage_variant_completeness = max(completeness)), by = .(MGG_K_locus, wgrr95)]

# Merge with K-locus summary for proper sorting
phage_variant_completeness <- merge(phage_variant_completeness, k_locus_summary, by = "MGG_K_locus")

p.completeness.pv95 <- ggplot(phage_variant_completeness, aes(x = reorder(MGG_K_locus, -total_SCs), y = phage_variant_completeness)) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "#D73027") +  # Red points for clarity
  geom_boxplot(outlier.shape = NA, fill = "#A6CEE3", color = "black", alpha = 0.6) + 
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of Phage Variant Completeness (WGRR95) Across K-loci",
    x = "K-locus",
    y = "Phage Variant Completeness (%)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

p.completeness.all.outfile <- "suppl-completeness-all.png"
p.completeness.p95.outfile <- "suppl-completeness-p95.png"

ggsave(p.completeness.all.outfile, p.completeness.all, width = 15, height = 7)
ggsave(p.completeness.p95.outfile, p.completeness.pv95, width = 15, height = 7)
