### load libraries
set.seed(3)
suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(ggsignif))
suppressWarnings(library(patchwork))


### load config
root <- rprojroot::find_rstudio_root_file()
config.filename <- file.path(root, "config", "config.yaml")
cfg <- yaml::read_yaml(config.filename)

# now set wd to the location of this script
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

phrogs.prob <- cfg$params$phrogs_prob
phrogs.qcov.fig5 <- cfg$params$phrogs_qcov_fig5
phrogs.tcov.fig5 <- cfg$params$phrogs_tcov_fig5


### LOAD DATA ###
# GWAS data
gwas.data.path <- file.path(cfg$paths$bogna$main, cfg$paths$bogna$raw_db_input_rel)
prophage.metadata.path = file.path(gwas.data.path,'/prophages_metadata.tsv')
pcs2proteins.dir = paste0(gwas.data.path, "/pcs2proteins.tsv")
raw.annots.dir = paste0(gwas.data.path, "/raw_hhsuite.tsv")

gwas.suppl.path <- file.path(cfg$paths$bogna$main, cfg$paths$bogna$supplement_data)
supl.table.st5.dir = paste0(gwas.suppl.path, "/SupplementaryTable_S5.xlsx")
supl.table.st6.dir = paste0(gwas.suppl.path, "/SupplementaryTable_S6.xlsx")


tableST5 = readxl::read_xlsx(supl.table.st5.dir)
tableST6 = readxl::read_xlsx(supl.table.st6.dir)
tested.depos.raw = tableST5 %>%
  rename(proteinID_wetlab = proteinID) %>%
  distinct(proteinID_wetlab, active) %>%
  inner_join(tableST6 %>% distinct(proteinID_wetlab, prophageID), by = "proteinID_wetlab") %>%
  mutate(active = if_else(active == "true", TRUE, if_else(active == "false", FALSE, NA))) 

# Original table with prophageId and proteinID and activity of tested depos used previously
# Note the order of proteins in the table is different from the supplementary one.
# gwas.depos.path <- file.path(cfg$paths$bogna$main, cfg$paths$mgg$db_input_rel)
# depos.prophages.dir = paste0(gwas.depos.path, '/depos_prophages.tsv')
# tested.depos.raw = fread(depos.prophages.dir) 
# bogna's local path: data.table::fread("/Users/bsmug/MGG Dropbox/Bogna Smug/Projects/2025_Phage_EcoEvo/Klebsiella_data/data/raw/depos_prophages.tsv") 


type.colors = c("expressed & active" = "darkred", "not expressed \n or not active" = "deepskyblue")
prophage_data = data.table::fread(prophage.metadata.path)


pcs2proteins = read.table(pcs2proteins.dir, header = TRUE) %>%
  distinct(PC, proteinID) %>%
  data.table::as.data.table()
pcs2proteins[, c("prophageID", "protID") := data.table::tstrsplit(proteinID, "_PROTEIN_", fixed = TRUE)]
pcs2proteins = pcs2proteins %>%
  mutate(protID = as.integer(protID))

annot_table = read.table(raw.annots.dir, 
                         sep = "\t", 
                         header = TRUE, 
                         quote = "",         # turn off quoting
                         comment.char = "",  # turn off comment interpretation
                         fill = TRUE) %>% 
  rename(PC = query) %>%
  filter(prob >= phrogs.prob)


phrog_annot_table = annot_table %>% 
  filter(db == "PHROGS" & qcov >= phrogs.qcov.fig5 & tcov >= phrogs.tcov.fig5) %>%
  rename(PHROG.category = category,
         PHROG.annotation = annot) %>%
  distinct(PC, PHROG.category, PHROG.annotation) %>%
  data.table::as.data.table() 


prophage.proteins.with.phrogs =  pcs2proteins %>% 
  left_join(phrog_annot_table, by = "PC", relationship = "many-to-many") %>%
  distinct(proteinID, prophageID, protID, PC, PHROG.annotation, PHROG.category) %>%
  mutate(PHROG.annotation = ifelse(is.na(PHROG.annotation), "unknown function", PHROG.annotation),
         PHROG.category = ifelse(is.na(PHROG.category), "unknown function", PHROG.category)) %>%
  data.table::as.data.table()



tested.depos = tested.depos.raw %>% 
  select(proteinID_wetlab, prophageID, active) %>%
  left_join(prophage_data, by = "prophageID") %>% 
  mutate(type = if_else(active,  "expressed & active", "not expressed \n or not active")) 

# Are there any signs phages with inactive depos are cryptic?
# Completeness
p2a =ggplot(data = tested.depos) + 
  geom_histogram(aes(x=completeness, fill = type), binwidth = 5,  boundary = 0, width = 0.5, col = "black") +
  theme_bw() + 
  xlim(c(0,100)) +
  ylab("Num. tested depolymerases") +
  xlab("Completeness") +
  scale_fill_manual(values = type.colors)
p2a
 tbl <- table(tested.depos$type, tested.depos$completeness)
fisher.test(tbl) 
xlsx::write.xlsx(x=tested.depos %>% select(proteinID_wetlab, completeness, type), file = "Data_FigS7.xlsx", sheetName = "A")

# Note geom_signif doesn't work with Fishers test when data has this format, so we will manually input the pre-computed p-values
p_to_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "NS.")
         )
  )
}



# NUMBER OF BACTERIAL ST
#Almost all prophages from which we isolated our 50 tested depolymerases can be found in 1 - 3 bacterial strains, regardless of the depolymerase activity
tested.phages.freq = tested.depos %>%
  distinct(proteinID_wetlab, type, wgrr95) %>%
  left_join(prophage_data %>% distinct(wgrr95, ST), by = "wgrr95", relationship = "many-to-many") %>%
  group_by(proteinID_wetlab, type) %>%
  summarise(num.ST = n_distinct(ST)) 

tbl <- table(tested.phages.freq$type, tested.phages.freq$num.ST)
ft1=fisher.test(tbl, alternative = "less") 
wt1=wilcox.test(num.ST ~ type, data = tested.phages.freq, alternative = "less")


p2=ggplot(tested.phages.freq %>% mutate(num.ST = as.integer(num.ST)), aes(x=type, y = num.ST, col = type)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  geom_jitter(height = 0, alpha = 0.5) +
  theme_bw() +
  #geom_signif(comparisons = list(c( "expressed & active" , "not expressed \n or not active")), 
  #            test = "wilcox.test", map_signif_level = TRUE, col = "black", step_increase = 0.1) +
  geom_signif(
    xmin = "expressed & active", xmax = "not expressed \n or not active",
    annotations = p_to_star(ft1$p.value),
    y_position = 8,
    color = "black"
  ) +
  guides(fill = "none") +
  xlab("") +
  ylab("Num. STs in which prophage cluster PV95 is detected in") +
  scale_y_continuous(breaks = scales::breaks_pretty(), limits = c(0,10)) +
  scale_color_manual(values = type.colors)
p2
xlsx::write.xlsx(x=tested.phages.freq %>% as.data.table(), file = "Data_FigS7.xlsx", sheetName = "B", append = TRUE)

## PHAGE CLUSTER SIZE
#There is no significant difference in phage cluster size for phages with active and inactive depos
wgrr95.clusters = prophage_data %>% 
  group_by(wgrr95) %>%
  summarise(num.prophages.in.cluster = n_distinct(prophageID),
            all.complete = all(completeness > 99.99))

cluster.size.data = tested.depos %>%
  distinct(proteinID_wetlab, prophageID, type) %>%
  left_join(prophage_data %>% distinct(prophageID, wgrr95)) %>%
  left_join(wgrr95.clusters)  

# Both Wilcox and fisher tests show no signficance p > 0.2, on plot we show wilcox but it is NS anyway
tbl <- table(cluster.size.data$type, cluster.size.data$num.prophages.in.cluster)
ft2=fisher.test(tbl, alternative = "less") 
wt2=wilcox.test(num.prophages.in.cluster ~ type, data = cluster.size.data, alternative = "less")

p3=ggplot(cluster.size.data, aes(x = type, y = num.prophages.in.cluster, col = type)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) + 
  geom_jitter(height=0) + 
  theme_bw() +
  #geom_signif(comparisons = list(c( "expressed & active" , "not expressed \n or not active")), 
  #            test = "wilcox.test", 
  #            map_signif_level = TRUE, 
  #            col = "black", 
  #            step_increase = 0.1) +
  geom_signif(
    xmin = "expressed & active", xmax = "not expressed \n or not active",
    annotations = p_to_star(ft2$p.value),
    y_position = 45,
    color = "black"
  ) +
  guides(fill = "none") +
  xlab("") +
  ylab("Num. prophages per prophage cluster PV95") +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  scale_color_manual(values = type.colors)
p3
xlsx::write.xlsx(x=cluster.size.data %>% 
                   select(proteinID_wetlab, type, prophageID, wgrr95, num.prophages.in.cluster) %>% as.data.frame(), 
                 file = "Data_FigS7.xlsx", sheetName = "C", append = TRUE)


## NUMBER OF TRANSPOSASES
transp.stats = tested.depos %>% 
  distinct(proteinID_wetlab, prophageID, wgrr95, type) %>%
  inner_join(prophage.proteins.with.phrogs, by = "prophageID", relationship = "many-to-many") %>%
  group_by(proteinID_wetlab, type, prophageID) %>%
  summarise(
    num.transposase = n_distinct(proteinID[PHROG.annotation == "transposase"])
  ) %>%
  mutate(num.transposase = as.integer(num.transposase))

# Both Wilcox and fisher tests show no signficance p > 0.2, on plot we show wilcox but it is NS anyway
tbl <- table(transp.stats$type, transp.stats$num.transposase)
ft3=fisher.test(tbl, alternative = "less") 
wt3=wilcox.test(num.transposase ~ type, data = transp.stats, alternative = "less")


p4=ggplot(transp.stats, aes(x = type, y = num.transposase, color = type)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  geom_jitter(height=0) + 
  theme_bw() +
  geom_signif(comparisons = list(c( "expressed & active" , "not expressed \n or not active")), 
              test = "wilcox.test", map_signif_level = TRUE, col = "black", step_increase = 0.1) +
  guides(fill = "none") +
  xlab("") +
  ylab("Num. transposases per prophage genome") +
  scale_y_continuous(breaks = seq(0:max(transp.stats$num.transposase) +1)) +
  scale_color_manual(values = type.colors)

p4
xlsx::write.xlsx(x=transp.stats %>% 
                   select(proteinID_wetlab, type, prophageID, num.transposase) %>%
                   as.data.frame(), 
                 file = "Data_FigS7.xlsx", sheetName = "D", append = TRUE)
#Glue plots
p2n=p2   + guides(color = "none")
p3n=p3   + guides(color = "none")
p4n=p4   + guides(color = "none")
p = (p2a | p2n) / (p3n | p4n) + plot_annotation(tag_levels = 'A') + plot_layout(guides = 'collect')
p
ggsave(filename =  'Activity_vs_cryptic.png', p, width = 25, height=25, units = "cm")
ggsave(filename =  'Activity_vs_cryptic.pdf', p, width = 25, height=25, units = "cm")


