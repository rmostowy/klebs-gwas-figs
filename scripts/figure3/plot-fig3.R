### load libraries
suppressWarnings(library(yaml))
suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))
suppressWarnings(library(ape))
suppressWarnings(library(ggtree))
suppressWarnings(library(ggtreeExtra))
suppressWarnings(library(ggplot2))
suppressWarnings(library(ggrepel))
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(forcats))
suppressWarnings(library(ggstar))
suppressWarnings(library(patchwork))

### load config
root <- rprojroot::find_rstudio_root_file()
config.filename <- file.path(root, "config", "config.yaml")
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

pyseer.filename <- file.path(data.path, cfg$data$gwas_hits_functs)
pyseer <- fread(pyseer.filename)
pyseer[ECOD1 == "", ECOD1 := "unknown"]
pyseer[ECOD2 == "", ECOD2 := "unknown"]

### LOAD PARMS ###

SC.CUTOFF.GWAS <- cfg$params$sc_cutoff_gwas


# Establish the GWAS-relevant K-loci
k_locus_diversity <- bacteria.metadata[, .(total_SCs = uniqueN(MGG_SC)), by = MGG_K_locus]
k_locus_diversity.gwas <- k_locus_diversity[total_SCs >= SC.CUTOFF.GWAS]
k_locus_diversity.gwas <- k_locus_diversity.gwas[MGG_K_locus != "UNK"]
k_locus_diversity.gwas[, k.index := as.numeric(gsub("KL", "", MGG_K_locus))]
setorder(k_locus_diversity.gwas, k.index)
Nk <- nrow(k_locus_diversity.gwas)

# Process GWAS results
PYSEER.CUTOFF <- -log10(0.05)
pyseer.results.lasso <- pyseer[minus_log10_pvalue_corr >= PYSEER.CUTOFF & beta > 0 & mode == "lasso"]
pyseer.results.lasso[, score := MCC*F1_score]

BEST.RESULT.PREC.THR <- cfg$params$best_results_prec_thr

# Algorithm for best result in Figure 3:
gwas.results.per.k.list <- lapply(1:Nk, function(index){
  # cat(index,"\n")
  # load all hits for the k locus
  this.k <- k_locus_diversity.gwas$MGG_K_locus[index]
  out <- data.table(k.locus = this.k, clustering = NA, PC = NA, mcc = NA, f1 = NA, precision = NA, recall = NA,
                    funct = NA, ecod1 = NA, ecod2 = NA)
  this.k.pyseer <- pyseer.results.lasso[locus == this.k]
  # filter by precision of BEST.RESULT.PREC.THR
  # this.k.pyseer.high.prec <- this.k.pyseer[precision_CI5_LOWER >= BEST.RESULT.PREC.THR]
  this.k.pyseer.high.prec <- this.k.pyseer[precision >= BEST.RESULT.PREC.THR & F1_score >= 0.5 & MCC >= 0.5]
  if(nrow(this.k.pyseer.high.prec) > 0){
    # of the remaining ones, pick the one with the highest F1*MCC score
    this.k.pyseer.best <- this.k.pyseer.high.prec[which.max(score)]
    out <- this.k.pyseer.best %>% select(k.locus = locus, clustering = version, PC, 
                                         mcc = MCC, f1 = F1_score, precision, recall,
                                         funct = PHROGS1, ecod1 = ECOD1, ecod2 = ECOD2)  
  }
  return(out)
})
gwas.results.per.k <- rbindlist(gwas.results.per.k.list)
gwas.results.per.k[is.na(ecod1), ecod1 := "unknown"]
gwas.results.per.k[is.na(ecod2), ecod2 := "unknown"]

######################
### Plot Figure 3A ###
######################

# Define custom shape mapping for clustering
# color_palette <- RColorBrewer::brewer.pal(9, "Set1")
col.palette <- brewer.pal(8, "Dark2")
col.palette2 <- brewer.pal(9, "Set1")
green.col <- col.palette[5]
yellow.col <- col.palette[6]
purple.col <- col.palette[3]
magenta.col <- col.palette[4]
turqoise.col <- col.palette[1]
brown.col <- col.palette[2]
beige.col <- col.palette[7]
pink.col <- col.palette2[8]
grey.col <- col.palette[8]

# Prepare the data
gwas.results.per.k.plot <- gwas.results.per.k %>%
  # filter(!is.na(k.locus), !is.na(recall)) %>%  # Remove missing values
  mutate(ecod1_adjusted = ifelse(grepl("Pectin lyase-like", ecod2), ecod2, ecod1))

pl.domain <- pyseer.results.lasso$ECOD1[grep("pectin", pyseer.results.lasso$ECOD1, ignore.case = T)[1]]
sgnh.domain <- pyseer.results.lasso$ECOD1[grep("sgnh", pyseer.results.lasso$ECOD1, ignore.case = T)[1]]
intra.domain <- pyseer.results.lasso$ECOD1[grep("intramolecular", pyseer.results.lasso$ECOD1, ignore.case = T)[1]]

nterminal.domain <- pyseer.results.lasso$ECOD1[grep("N-terminal", pyseer.results.lasso$ECOD1, ignore.case = T)[1]]
sixbladed.domain <- pyseer.results.lasso$ECOD1[grep("6-bladed", pyseer.results.lasso$ECOD1, ignore.case = T)[1]]
ploop.domain <- pyseer.results.lasso$ECOD1[grep("P-loop", pyseer.results.lasso$ECOD1, ignore.case = T)[1]]
nad.domain <- pyseer.results.lasso$ECOD1[grep("NAD/FAD", pyseer.results.lasso$ECOD1, ignore.case = T)[1]]
galactose.domain <- pyseer.results.lasso$ECOD1[grep("Galactose", pyseer.results.lasso$ECOD1, ignore.case = T)[1]]

unknown.domain <- "unknown"

# Redefine color mapping dynamically
ecod_all_colors <- setNames(
  c(green.col, yellow.col,  purple.col,   brown.col,        beige.col,              turqoise.col,              magenta.col,  pink.col,   grey.col),
  c(pl.domain, sgnh.domain, intra.domain, nterminal.domain, sixbladed.domain, galactose.domain, ploop.domain, nad.domain, unknown.domain)
)

gwas.results.per.k.plot$k.locus <- factor(gwas.results.per.k.plot$k.locus, levels = k_locus_diversity.gwas$MGG_K_locus)
# gwas.results.per.k.plot$clustering <- factor(gwas.results.per.k.plot$clustering, exclude = NA)
# Define shape mapping for clustering
# shape_mapping <- c("PCI80C80" = 22, "PCI50C50" = 23, "PCI80C50" = 24, "PCI50C80" = 0, "PCI00C80" = 5, "PCI00C50" = 2) 
# Define star shapes for clustering categories
star_shapes <- c("PCI80C80" = 6, "PCI50C50" = 14, "PCI80C50" = 9, 
                 "PCI50C80" = 13, "PCI00C80" = 11, "PCI00C50" = 15)

plot.gwas <- ggplot(gwas.results.per.k.plot, aes(x = k.locus, y = recall, color = ecod1_adjusted, shape = clustering)) +
  geom_star(aes(fill = ecod1_adjusted, starshape = clustering), size = 8, alpha = 0.8, color = "black") +   # Ensures black outline on all shapes
  scale_fill_manual(values = ecod_all_colors, na.translate = FALSE) +
  scale_color_manual(values = ecod_all_colors, na.translate = FALSE) +  
  scale_starshape_manual(values = star_shapes, na.translate = FALSE) +  
  theme_minimal(base_size = 14) +  
  labs(x = "K Locus",  
       y = "Recall",  
       starshape = "Clustering Level",   # Custom legend title for clustering
       fill = "ECOD Domains") +  # Custom legend title for ECOD1 annotations
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        legend.key = element_blank()) +  # Removes black rectangles from legend
  scale_y_continuous(limits = c(0, 1)) + 
  guides(
    starshape = guide_legend(
      override.aes = list(size = 5, size = 4, color = "black")
    ),
    fill = guide_legend(
      override.aes = list(starshape = 15, size = 5, color = "black")
    )  # Changes ECOD fill legend to a circle (shape 21)
  )
ggsave("Figure3A.pdf", plot.gwas, width = 16, height = 5)


######################
### Plot Figure 3B ###
######################

# By threshold
prec.thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)
ecods.vs.precision.list <- lapply(1:length(prec.thresholds), function(index){
  this.prec.thr <- prec.thresholds[index]
  this.prec.thr.ecods.list <- lapply(1:Nk, function(k.index){
    out.ecod <- NULL
    this.k <- k_locus_diversity.gwas$MGG_K_locus[k.index]
    this.k.pyseer <- pyseer.results.lasso[locus == this.k]
    this.k.pyseer.high.prec <- this.k.pyseer[precision >= this.prec.thr & F1_score >= 0.5 & MCC >= 0.5]
    if(nrow(this.k.pyseer.high.prec) > 0){
      this.k.pyseer.best.clustering <- this.k.pyseer.high.prec[which.max(score)]$version
      this.k.pyseer.best <- this.k.pyseer.high.prec[version == this.k.pyseer.best.clustering]
      
      all.ecod.hits.list <- lapply(1:nrow(this.k.pyseer.best), function(i){
        out <- "unknown"
        ecod.prediction1 <- this.k.pyseer.best[i]$ECOD1
        ecod.prediction2 <- this.k.pyseer.best[i]$ECOD2
        ecod.predictions <- c(ecod.prediction1, ecod.prediction2)
        ecod.predictions.known <- ecod.predictions[ecod.predictions != "unknown"]
        if(length(ecod.predictions.known) > 0) out <- unique(ecod.predictions.known)
        out  
      })
      all.ecod.hits <- unlist(all.ecod.hits.list)
      out.ecod <- all.ecod.hits
    }
    out.ecod
  })
  this.prec.thr.ecods <- unlist(this.prec.thr.ecods.list)
  data.table(min.precision = this.prec.thr, ecod.domain = this.prec.thr.ecods)
})
ecods.vs.precision <- rbindlist(ecods.vs.precision.list)

#######################
# By score
score.thresholds <- c(0.1, 0.2, 0.3, 0.4, 0.5)
ecods.vs.score.list <- lapply(1:length(score.thresholds), function(index){
  this.score.thr <- score.thresholds[index]
  this.prec.thr.ecods.list <- lapply(1:Nk, function(k.index){
    out.ecod <- NULL
    this.k <- k_locus_diversity.gwas$MGG_K_locus[k.index]
    this.k.pyseer <- pyseer.results.lasso[locus == this.k]
    this.k.pyseer.high.prec <- this.k.pyseer[precision >= 0.8 & F1_score >= this.score.thr & MCC >= this.score.thr]
    if(nrow(this.k.pyseer.high.prec) > 0){
      this.k.pyseer.best.clustering <- this.k.pyseer.high.prec[which.max(score)]$version
      this.k.pyseer.best <- this.k.pyseer.high.prec[version == this.k.pyseer.best.clustering]
      
      all.ecod.hits.list <- lapply(1:nrow(this.k.pyseer.best), function(i){
        out <- "unknown"
        ecod.prediction1 <- this.k.pyseer.best[i]$ECOD1
        ecod.prediction2 <- this.k.pyseer.best[i]$ECOD2
        ecod.predictions <- c(ecod.prediction1, ecod.prediction2)
        ecod.predictions.known <- ecod.predictions[ecod.predictions != "unknown"]
        if(length(ecod.predictions.known) > 0) out <- unique(ecod.predictions.known)
        out  
      })
      all.ecod.hits <- unlist(all.ecod.hits.list)
      out.ecod <- all.ecod.hits
    }
    out.ecod
  })
  this.prec.thr.ecods <- unlist(this.prec.thr.ecods.list)
  data.table(min.score = this.score.thr, ecod.domain = this.prec.thr.ecods)
})
ecods.vs.score <- rbindlist(ecods.vs.score.list)

#######################
# Aggregate counts for both datasets
ecods.vs.precision.count <- ecods.vs.precision %>% count(min.precision, ecod.domain)
ecods.vs.score.count <- ecods.vs.score %>% count(min.score, ecod.domain)

# Rename the precision column to a generic name
setnames(ecods.vs.precision.count, "min.precision", "Threshold")
setnames(ecods.vs.score.count, "min.score", "Threshold")

# Add a new column to distinguish between Precision and Score
ecods.vs.precision.count$Threshold_Type <- "Precision"
ecods.vs.score.count$Threshold_Type <- "Score"

# Combine both datasets into one
ecods.combined <- rbind(ecods.vs.precision.count, ecods.vs.score.count)

# Ensure Threshold is treated as a factor for correct faceting order
ecods.combined$Threshold <- factor(ecods.combined$Threshold, levels = sort(unique(ecods.combined$Threshold)))

# Create custom facet labels to control ordering
ecods.combined$Facet_Label <- factor(
  paste(ecods.combined$Threshold_Type, ecods.combined$Threshold, sep = " ≥ "),
  levels = c(
    paste("Precision", sort(unique(ecods.vs.precision$min.precision)), sep = " ≥ "),
    paste("Score", sort(unique(ecods.vs.score$min.score)), sep = " ≥ ")
  )
)

# Order ECOD domains based on total counts across both datasets
domain_order <- ecods.combined %>%
  group_by(ecod.domain) %>%
  summarise(total_count = sum(n)) %>%
  arrange(desc(total_count)) %>%
  pull(ecod.domain)

# Apply this ordering to the dataset
ecods.combined$ecod.domain <- factor(ecods.combined$ecod.domain, levels = rev(domain_order))

COUNT.THR <- 3
ecods.total <- ecods.combined %>%
  group_by(ecod.domain) %>%
  summarise(total_n = sum(n), .groups = "drop")
setDT(ecods.total)
filtered.domains <- as.character(ecods.total[total_n >= COUNT.THR]$ecod.domain)

ecods.combined <- ecods.combined[ecod.domain %in% filtered.domains]

# ✅ Ensure all ECOD domains have a color (default black if missing)
ecod_all_colors_full <- ecod_all_colors  # Copy the original color dictionary
missing_domains <- setdiff(unique(ecods.combined$ecod.domain), names(ecod_all_colors_full))
for (domain in missing_domains) {
  ecod_all_colors_full[domain] <- "#000000"
}


#######################
# plot data

# ✅ Generate the improved 2-row, 5-column faceted plot
plot.ecods.combined <- ggplot(ecods.combined, 
                              aes(y = fct_reorder(ecod.domain, n), x = n, fill = ecod.domain)) +
  geom_bar(stat = "identity") +
  
  # ✅ Facet using the new ordered `Facet_Label` to enforce correct row-column structure
  facet_wrap(~ Facet_Label, nrow = 2) +
  
  theme_minimal(base_size = 14) +
  
  # ✅ Custom fill colors, disable legend
  scale_fill_manual(values = ecod_all_colors_full, na.translate = FALSE, guide = "none") +
  
  labs(
    x = "Count",  
    y = ""
  ) +
  
  theme(
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 10, face = "bold"),  
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

# Display the final combined plot
print(plot.ecods.combined)
ggsave("Figure3B.pdf", plot.ecods.combined, width = 16, height = 5)
