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
pyseer.filename <- file.path(data.path, "pyseer_hits_all.csv")
pyseer <- fread(pyseer.filename)

funct.table.filename <- file.path(data.path, "clusters_functions_best_all.csv")
funct.table <- fread(funct.table.filename)

ecod.domains <- fread(file.path(root, 'scripts', 'data', cfg$data$ecod_domains))
setkey(ecod.domains, full_id)

### LOAD PARMS ###

phrogs.prob <- cfg$params$phrogs_prob
phrogs.qcov <- cfg$params$phrogs_qcov
phrogs.scov <- cfg$params$phrogs_scov

ecod.prob <- cfg$params$ecod_prob
ecod.qcov <- cfg$params$ecod_qcov
ecod.scov <- cfg$params$ecod_scov


### PROCESS DATA ###

# Filter functional hits
funct.table.filtered.phrogs <- funct.table[db == "PHROGS" & prob >= phrogs.prob & qcov >= phrogs.qcov & tcov >= phrogs.scov & annot != "unknown function"]
funct.table.filtered.ecod <- funct.table[db == "ECOD" & prob >= ecod.prob & qcov >= ecod.qcov & tcov >= ecod.scov]

funct.table.filtered.ecod$t.id <- ecod.domains[.(funct.table.filtered.ecod$target)]$t_id
funct.table.filtered.ecod$t.name <- ecod.domains[.(funct.table.filtered.ecod$target)]$t_name

# Pick best PHROG hits
best_per_function <- funct.table.filtered.phrogs[
  ,
  .SD[which.max(bits)],
  by = .(version, PC, annot)
]
best_per_function <- best_per_function[
  order(-bits)
][
  ,
  rank := seq_len(.N),
  by = .(version, PC)
]
top2_functions <- best_per_function[rank <= 2, .(version, PC, rank, annot)]
phrogs.top.two <- dcast(
  top2_functions,
  version + PC ~ rank,
  value.var = "annot"
)
setnames(phrogs.top.two, c("1", "2"), c("PHROGS1", "PHROGS2"))


# Pick best ECOD hits
best_per_domain <- funct.table.filtered.ecod[
  ,
  .SD[which.max(bits)],
  by = .(version, PC, t.name)
]

best_per_domain <- best_per_domain[
  order(-bits)
][
  ,
  rank := seq_len(.N),
  by = .(version, PC)
]

top2_domains <- best_per_domain[rank <= 2, .(version, PC, rank, t.name)]

ecods.top.two <- dcast(
  top2_domains,
  version + PC ~ rank,
  value.var = "t.name"
)
setnames(ecods.top.two, c("1", "2"), c("ECOD1", "ECOD2"))

# Merge functional hits with the pyseer table
pyseer.funct <- merge(pyseer, phrogs.top.two, by = c("version", "PC"), all.x = TRUE)
pyseer.funct <- merge(pyseer.funct, ecods.top.two, by = c("version", "PC"), all.x = TRUE)
pyseer.funct.outfile <- file.path(data.path, cfg$data$gwas_hits_functs)
fwrite(pyseer.funct, file = pyseer.funct.outfile)


# Filter GWAS hits
PYSEER.CUTOFF <- -log10(0.05)
pyseer.default <- pyseer.funct[version == "PCI50C50" & mode == "lasso"]
pyseer.default.hits <- pyseer.default[minus_log10_pvalue_corr >= PYSEER.CUTOFF & beta > 0]

pyseer.results <- pyseer.default.hits %>% select(locus, PC, beta, pc.freq = PC_freq,
                                                 f1 = F1_score, f1.sc = F1_score_SC, 
                                                 mcc = MCC, mcc.sc = MCC_SC,
                                                 precision, precision.lower = precision_CI5_LOWER, precision.upper = precision_CI95_UPPER,
                                                 recall, recall.lower = recall_CI5_LOWER, recall.upper = recall_CI95_UPPER,
                                                 precision.sc = precision_SC, precision.sc.lower = precision_SC_CI5_LOWER, precision.sc.upper = precision_SC_CI95_UPPER,
                                                 recall.sc = recall_SC, recall.sc.lower = recall_SC_CI5_LOWER, recall.sc.upper = recall_SC_CI95_UPPER,
                                                 funct = PHROGS1, ecod1 = ECOD1, ecod2 = ECOD2
)
pyseer.results[is.na(funct), funct := "unknown"]
pyseer.results[funct == "unknown function", funct := "unknown"]
pyseer.results[, label := sprintf("%s (%s)", locus, PC)]

pyseer.results.filtered <- pyseer.results[f1 > 0.5 & mcc > 0.5]

# Plot the results

pyseer.results.filtered.plot <- pyseer.results.filtered
pyseer.results.filtered.plot[, label2 := label]
pyseer.results.filtered.plot[f1 < 0.7 & mcc < 0.7]$label2 <- ""

gwas.results <- ggplot(pyseer.results.filtered.plot, 
                       aes(x = f1, y = mcc, label = label2, color = funct)) +
  geom_point(alpha = .9, size = 3) +
  geom_text_repel(size = 5, box.padding = 0.5, max.overlaps = 20, colour = "black") +
  scale_color_brewer(palette = "Paired") +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = "black") +
  theme_bw(base_size = 14) +
  labs(
    x = "F1 Score",
    y = "Matthews Correlation Coefficient (MCC)",
    color = "Function (PHROGs)"
  ) +
  
  scale_x_continuous(limits = c(0.4, 1)) +
  scale_y_continuous(limits = c(0.5, 1))

gwas.results.pr <- ggplot(pyseer.results.filtered.plot, aes(x = precision, y = recall, label = label, colour = funct)) + 
  geom_smooth(method = "lm", colour = "black", fullrange = TRUE) +  # Fit a linear model
  geom_point(alpha = 0.9, size = 3) +          # Scatter points
  geom_text_repel(size = 3, box.padding = 0.5, max.overlaps = 20, colour = "black") +
  scale_color_brewer(palette = "Paired") +
  # geom_text_repel(aes(label = PC), size = 4, box.padding = 0.5, max.overlaps = 20) +
  scale_x_continuous(limits = c(0.25, 1)) + 
  scale_y_continuous(limits = c(0.25, 1)) + 
  theme_bw(base_size = 14) + 
  guides(color = "none") +
  labs(
    x = "Precision",
    y = "Recall",
    color = "Function (PHROGs)"
  )

figure.ab <- gwas.results + gwas.results.pr + plot_annotation(tag_levels = 'A')
# ggsave("Figure2A_B.png", figure.ab, width = 12, height = 5)
ggsave("Figure2A_B.pdf", figure.ab, width = 18, height = 8)
