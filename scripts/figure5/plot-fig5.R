### load libraries
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(patchwork))
suppressWarnings(library(UpSetR))


### load config
root <- rprojroot::find_rstudio_root_file()
config.filename <- file.path(root, "config", "config.yaml")
cfg <- yaml::read_yaml(config.filename)

# now set wd to the location of this script
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

####################################################################
############################### MAIN ###############################
####################################################################
### LOAD PARAMS ###
phrogs.prob <- cfg$params$phrogs_prob
phrogs.qcov.fig5 <- cfg$params$phrogs_qcov_fig5
phrogs.tcov.fig5 <- cfg$params$phrogs_tcov_fig5

ecod.prob <- cfg$params$ecod_prob
ecod.tcov.fig5 <- cfg$params$ecod_tcov_fig5
min.completeness = cfg$params$min_completeness

### DEFINE ADDITIONAL PARAMS RELATED TO RESULTS FILTERING AND VISUALISATIION ###
# Parameters (whcih clusters and ECODs to plot)
MIN.CLUST.SIZE = 10
MIN.NUM.SC = 3
TEXT_SIZE = 14
# Define colours
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
blue.col <- col.palette2[2]

ECOD_labeller = function(ecods) {
  from = c("SGNH hydrolase", "other",
           "Multiple","Pectin lyase-like",                                
           "Phage tail fiber protein trimerization domain",
           "gp11/gp12 receptor-binding domain",                 
           "Putative tailspike protein Orf210 N-terminal domain",
           "Intramolecular chaperone domain in virus tail spike protein")
  to = c("SGNH hydrolase", "other",
         "Multiple","Pectin lyase-like",                                
         "Phage tail fiber protein\ntrimerization domain",
         "gp11/gp12 receptor-binding domain",                 
         "Putative tailspike protein\nOrf210 N-terminal domain",
         "Intramolecular chaperone domain\nin virus tail spike protein")
  return(plyr::mapvalues(ecods, from, to, warn_missing = TRUE))
} 

# Define color mapping for selected ECODs (after manual verification of the most frequent domains)
custom_colors <- setNames(
  c(green.col, yellow.col, purple.col, 
    brown.col, 
    "lightcoral",
    magenta.col,
    "lightblue",
    "darkviolet",
    "gray80", "gray30", "gray90"),
  c("Pectin lyase-like", "SGNH hydrolase",  "Intramolecular chaperone domain in virus tail spike protein",
    "Putative tailspike protein Orf210 N-terminal domain",
    "Phage tail fiber protein trimerization domain",
    "gp11/gp12 receptor-binding domain", 
    "Lysozyme-like" ,
    "Uncharacterized protein CV0426",
    "None of above", "Multiple", "No RBP detected"))

names(custom_colors) = ECOD_labeller(names(custom_colors))

### LOAD DATA ###
# GWAS data
gwas.data.path <- file.path(cfg$paths$bogna$main, cfg$paths$mgg$db_input_rel)

gwas.metadata.filename <- file.path(gwas.data.path,"bacteria.tsv")
gwas.metadata <- fread(gwas.metadata.filename)

gwas.prophage.filename <- file.path(gwas.data.path, "prophages.tsv")
gwas.prophages <- fread(gwas.prophage.filename) %>% 
  distinct(prophageID)

# RAW PROPHAGE DATA
raw.prophage.data.path <- file.path(cfg$paths$bogna$main, cfg$paths$bogna$raw_db_input_rel)

raw.annots.filename <- file.path(raw.prophage.data.path, "raw_hhsuite.tsv")
annot_table = read.table(raw.annots.filename, 
                         sep = "\t", 
                         header = TRUE, 
                         quote = "",         # turn off quoting
                         comment.char = "",  # turn off comment interpretation
                         fill = TRUE)  %>% 
  rename(PC = query) 

prophages.filename <- file.path(raw.prophage.data.path, "prophages_metadata.tsv")
prophages = read.table(prophages.filename, header = TRUE, sep = "\t") 

pcs2proteins.filename <- file.path(raw.prophage.data.path, "pcs2proteins.tsv")
pcs2proteins = read.table(pcs2proteins.filename, header = TRUE) %>%
  distinct(PC, proteinID) %>%
  as.data.table()

### PROCESS DATA ###
setkey(gwas.metadata, genomeID)
prophages_table = prophages %>%
  left_join(gwas.metadata %>% distinct(genomeID, MGG_SC), by = "genomeID")

quality_prophages_table = prophages_table %>% 
  filter(confidence == "high" & completeness > min.completeness) %>%
  # only use prophages used in GWWAS
  inner_join(gwas.prophages, by = "prophageID") %>%
  distinct(prophageID, MGG_K_locus)

pcs2proteins[, c("prophageID", "protID") := tstrsplit(proteinID, "_PROTEIN_", fixed = TRUE)]
pcs2proteins = pcs2proteins %>%
  mutate(protID = as.integer(protID))

# Extract most diverse K-loci from GWAS
k_locus_diversity <- gwas.metadata[, .(total_SCs = uniqueN(MGG_SC)), by = MGG_K_locus]
k_locus_diversity <- k_locus_diversity[MGG_K_locus != "UNK"]
top.k.loci <- k_locus_diversity[total_SCs >= 10]
top.k.loci[, locus.index := as.numeric(gsub("KL", "", MGG_K_locus))]
top.k.loci <- top.k.loci[order(locus.index)]
k.locus.order <- top.k.loci$MGG_K_locus

# Get, filter and clean ECOD annotations
ecod_annot_table = annot_table %>% 
  filter(db == "ECOD" & tcov >= ecod.tcov.fig5 & prob >= ecod.prob) %>%
  distinct(PC, annotation = name, target) %>%
  as.data.table()
ecod_annot_table[, c("ecod", "#uid", "ecod_domain_id") := tstrsplit(target, "_", fixed=TRUE)]
ecod_annot_table[, c("name", "number", "name2", "ecod.level.desc", "ecod.description") := tstrsplit(annotation, "|", fixed=TRUE)]
ecod_annot_table[, c("axh", "ecod.tf") := tstrsplit(ecod.level.desc, ", T: ", fixed=TRUE)]
ecod_annot_table[, c("ecod.t", "ecod.f") := tstrsplit(ecod.tf, ", F: ", fixed=TRUE)]
ecod_annot_table_f = ecod_annot_table %>% distinct(PC, axh, ecod.t, ecod.f)
ecod_annot_table_t = ecod_annot_table %>% distinct(PC, ecod.t)

# Get, filter and clean PHROG annotations
phrog_annot_table = annot_table %>% 
  filter(db == "PHROGS" & qcov >= phrogs.qcov.fig5 & tcov >= phrogs.tcov.fig5 & prob >= phrogs.prob) %>%
  rename(PHROG.category = category,
         PHROG.annotation = annot) %>%
  distinct(PC, PHROG.category, PHROG.annotation) %>%
  as.data.table() 


# Create a table of annotated prophage proteins 
# NOTE THAT PROTEINS WITH MULTIPLE ANNOTATIONS WILL APEAR HERE MULTIPLE TIMES
prophage.proteins.with.ecods =  pcs2proteins %>% 
  left_join(ecod_annot_table_t, by = "PC", relationship = "many-to-many") %>%
  distinct(proteinID, prophageID, protID, PC, ecod.t) %>%
  mutate(ecod.t = ifelse(is.na(ecod.t), "unknown function", ecod.t)) %>%
  as.data.table()

prophage.proteins.with.phrogs =  pcs2proteins %>% 
  left_join(phrog_annot_table, by = "PC", relationship = "many-to-many") %>%
  distinct(proteinID, prophageID, protID, PC, PHROG.annotation, PHROG.category) %>%
  mutate(PHROG.annotation = ifelse(is.na(PHROG.annotation), "unknown function", PHROG.annotation),
         PHROG.category = ifelse(is.na(PHROG.category), "unknown function", PHROG.category)) %>%
  as.data.table()

# CHECK IF THE FILTERING WAS SENSIBLE: HOW MANY FUNCTIONS WE GET PER PC
# how many functions per PC
pstats1 = phrog_annot_table %>% 
  filter( PHROG.annotation != "unknown function") %>%
  group_by(PC) %>%
  summarise(n.annot.per.pc = n_distinct(PHROG.annotation))
print(summary(pstats1$n.annot.per.pc))

# how many ECOD domains per PC
pstats2=prophage.proteins.with.ecods %>% 
  group_by(PC) %>%
  summarise(n.ecod.per.pc = n_distinct(ecod.t))
print(summary(pstats2$n.ecod.per.pc))



########################################################################
# 1. For each K locus within the top loci find all prophages in bacterial isolates with this K locus and find potential RBPs in these phages based on PHROG annotation
#rbp.grep.string <- "tail fiber protein|tail spike protein"
rbp.grep.string <- "tail fiber protein|tail spike protein|tail collar fiber protein|tail protein and host specificity|central tail fiber J"
# that also includes: # "lytic tail fiber protein","tail spike protein with colonic acid degradation activity"

all.rbp.associated.pcs.list <- lapply(1:nrow(top.k.loci), function(index){
  cat(index,"\n")
  this.k.locus = top.k.loci[index]$MGG_K_locus
  # 1. Find all prophages in bacterial isolates with this K locus
  this.k.locus.prophages = prophages_table %>% 
    filter(MGG_K_locus == this.k.locus) %>% 
    pull(prophageID)
  # 2. Find all proteins in those prophages (and their annptations)
  this.k.locus.prophage.functs = prophage.proteins.with.phrogs %>% 
    filter(prophageID %in% this.k.locus.prophages)
  # 3. Shortlist proteins with RBP hits by PHROGs
  is.rbp = grepl(rbp.grep.string, this.k.locus.prophage.functs$PHROG.annotation)
  this.k.locus.rbp.proteins = this.k.locus.prophage.functs[is.rbp]$proteinID %>% unique()
  # 4. Identify all annotation clusters which these RBPs belong to
  this.k.locus.rbp.clusters <- pcs2proteins %>% 
    filter(proteinID %in% this.k.locus.rbp.proteins) %>% 
    pull(PC) %>% 
    unique()
  return(this.k.locus.rbp.clusters)
})
all.rbp.associated.pcs.unique <- unlist(all.rbp.associated.pcs.list) %>% unique()

# Find all ECODs associated with our clusters
rbp.associated.ecods = prophage.proteins.with.ecods %>% 
  filter(PC %in% all.rbp.associated.pcs.unique & ecod.t != "unknown function") %>% 
  distinct(ecod.t) %>% 
  pull(ecod.t)

# For each of these ECOD find all clusters that may be described by this ECOD and calculate their diversity 
rbp.clust.data = data.table(cluster = character(0), size = integer(0), num.sc = integer(0), k.locus.simpson = numeric(0), ecod = character(0))
for (this.rbp.associated.ecod in rbp.associated.ecods) {
  # take all PCs with good hits to RBPs and with at least one protein in the PC annotated by this ECOD
  this.ecod.rbp.pcs = prophage.proteins.with.ecods %>% 
    filter(PC %in% all.rbp.associated.pcs.unique) %>%
    filter(ecod.t == this.rbp.associated.ecod) %>%
    pull(PC) %>%
    unique()
  no.clust <- length(this.ecod.rbp.pcs)
  # now for each of these clusters calculate its diversity of loci, size, num associated SC
  rbp.clust.data.list <- lapply(1:no.clust, function(cluster.index){
    this.cluster <- this.ecod.rbp.pcs[cluster.index]
    this.cluster.prophages <- pcs2proteins %>% filter(PC == this.cluster) %>% pull(prophageID) %>% unique()
    this.cluster.proteins <- pcs2proteins %>% filter(PC == this.cluster) %>% pull(proteinID) %>% unique()
    this.cluster.k.loci <- prophages_table %>% filter(prophageID %in% this.cluster.prophages) %>% pull(MGG_K_locus)
    this.cluster.k.loci <- this.cluster.k.loci[this.cluster.k.loci != "UNK"]
    this.cluster.k.loci.unique <- unique(this.cluster.k.loci)
    this.cluster.k.loci.unique.index <- as.numeric(gsub("KL", "", this.cluster.k.loci.unique))
    sorted.k.loci <- this.cluster.k.loci.unique[sort(this.cluster.k.loci.unique.index, index.return = T)$ix]
    
    data.table(cluster = this.cluster,
               size = length(this.cluster.proteins),
               num.sc = prophages_table %>% filter(prophageID %in% this.cluster.prophages) %>% distinct(MGG_SC) %>% nrow(),
               k.locus.simpson = vegan::diversity(table(this.cluster.k.loci), index = "simpson"),
               ecod = this.rbp.associated.ecod)
  })
  rbp.clust.data.this.ecod <- rbindlist(rbp.clust.data.list)
  rbp.clust.data = rbind(rbp.clust.data, rbp.clust.data.this.ecod)
}



# Plot data
# Only show top ECOD annotations sorting by the number of PCs where these ECODs were top annotations
rbp.clust.filtered = rbp.clust.data %>% 
  as.data.frame() %>%
  filter(size >= MIN.CLUST.SIZE & num.sc >= MIN.NUM.SC) %>% 
  group_by(ecod) %>% 
  mutate(num.pc = n_distinct(cluster)) %>%
  ungroup() %>%
  filter(num.pc > 1)

top.ecods = rbp.clust.filtered %>% 
  filter(ecod != "unknown function") %>% 
  group_by(ecod) %>%
  summarise(num.clusters = n_distinct(cluster)) %>%
  arrange(num.clusters) %>%
  distinct(ecod) %>% 
  pull(ecod) 


rbp.clust.filtered.to.plot = rbp.clust.filtered %>%
  mutate(ecod = ECOD_labeller(ecod),
         `# SC` = num.sc,
         ecod = factor(ecod, levels = ECOD_labeller(top.ecods)))

# K-locus diversity per ECOD domain
# SHowing only for most ppular ECODs, PCs which appeared in atleast 3 SCs and contain at least 20 proteins so that we can calculate the diversity index
plot.gen.spec <- ggplot(rbp.clust.filtered.to.plot, 
                        aes(x = ecod, y = k.locus.simpson, color = ecod)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(size = `# SC`), width = 0.2, height = 0,alpha = 0.7) + 
  scale_color_manual(values = ECOD_labeller(custom_colors), na.translate = FALSE) +  
  scale_size(range = c(4, 10), name = "Num. SC") +
  theme_bw(base_size = 14) + 
  labs(x = "", y = "K-locus diversity\n(Simpson index)",
       color = "ECOD Domain at RBP") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = TEXT_SIZE),  # Rotate X-axis labels
    axis.text.y = element_text(size = TEXT_SIZE),
    legend.title = element_text(size = TEXT_SIZE),
    legend.position = "right",
    legend.text = element_text(size = TEXT_SIZE)
  ) +
  guides(color = guide_legend(override.aes = list(size=6)))



#################################################
# Fig 5b, now look at all complete prophages regardless of what data we had in GWAS
# All prophage data, not only GWAS proteins

# protID are the integer numbers: 1,2,3,... etc denoting the order of proteins in phage genome
is.near.RBP = function(is.rbp, protID) {
  if (length(is.rbp) != length(protID)) {
    error("Provided vectors must be of the same length")
  }
  n = length(is.rbp)
  near.rbp = rep(FALSE, n)
  if (any(is.rbp)) {
    # find the protIDs of the proteins at RBP
    protIDRBP = protID[is.rbp]
    # Find protIDs of the proteins around RBP and make sure they are non negative
    prot.id.near.rbp = c(protIDRBP, protIDRBP + 1, protIDRBP + 2, protIDRBP - 1, protIDRBP - 2) %>% unique()  %>% sort()
    prot.id.near.rbp = prot.id.near.rbp[which(prot.id.near.rbp > 0)]
    idx.near.rbp = which(protID %in% prot.id.near.rbp)
    near.rbp[idx.near.rbp] = TRUE
  } 
  return(near.rbp)
}


Get.Phage.Type.Top.Ecods = function(phage_proteins, top.ecods, loc = "any", all.ecod.combinations = FALSE) {
    if (loc == "at rbp") {
      phage_proteins.with.type = phage_proteins %>%
        group_by(prophageID, MGG_K_locus, has.rbp) %>%
        summarise(phage.type = paste(sort(unique(ecod.t[ecod.t %in% top.ecods & is.rbp])), collapse = " & "),
               num.ecod.t = length(unique(ecod.t[ecod.t %in% top.ecods & is.rbp]))) %>%
        ungroup()
    } else if (loc == "around rbp") {
      phage_proteins.with.type = phage_proteins %>%
        group_by(prophageID, MGG_K_locus, has.rbp) %>%
        summarise(phage.type = paste(sort(unique(ecod.t[ecod.t %in% top.ecods & is.around.rbp])), collapse = " & "),
                  num.ecod.t = length(unique(ecod.t[ecod.t %in% top.ecods & is.around.rbp]))) %>%
        ungroup()
    } else if (loc == "any") {
      phage_proteins.with.type = phage_proteins %>%
        group_by(prophageID, MGG_K_locus, has.rbp) %>%
        summarise(phage.type = paste(sort(unique(ecod.t[ecod.t %in% top.ecods])), collapse = " & "),
               num.ecod.t = length(unique(ecod.t[ecod.t %in% top.ecods])))%>%
        ungroup()
    }
    if (!all.ecod.combinations) {
      ind = which(phage_proteins.with.type$num.ecod.t > 1)
      phage_proteins.with.type$phage.type[ind] = "Multiple"
    }
  ind.no.rbp = which(phage_proteins.with.type$phage.type == "" & !phage_proteins.with.type$has.rbp)
  phage_proteins.with.type$phage.type[ind.no.rbp] = "No RBP detected"
  ind = which(phage_proteins.with.type$phage.type == "" & phage_proteins.with.type$has.rbp)
  phage_proteins.with.type$phage.type[ind] = "None of above"
  return(phage_proteins.with.type)
} 

quality.phage.types = data.frame(prophageID = character(0), MGG_K_locus = character(0), phage.type = character(0))
for (this.k.locus in k.locus.order) {
  phages.this.k.locus = quality_prophages_table %>% 
    filter(MGG_K_locus == this.k.locus) %>%
    distinct(prophageID, MGG_K_locus)
  quality.phage.proteins.this.k.locus = prophage.proteins.with.ecods[prophageID %in% phages.this.k.locus$prophageID] %>%
    mutate(is.rbp = PC %in% all.rbp.associated.pcs.unique) %>%
    group_by(prophageID) %>%
    mutate(has.rbp = any(is.rbp)) %>%
    mutate(is.around.rbp = is.near.RBP(is.rbp, protID)) %>%
    ungroup() %>%
    select(prophageID, protID, ecod.t, is.rbp, is.around.rbp, has.rbp) %>%
    mutate(type = "other") %>%
    mutate(MGG_K_locus = this.k.locus)
  #setDT(quality.phage.proteins.this.k.locus)
  phage.types.this.k.locus.at.rbp = Get.Phage.Type.Top.Ecods(phage_proteins = quality.phage.proteins.this.k.locus, top.ecods = top.ecods, loc =  "at rbp", all.ecod.combinations = FALSE) %>%
    distinct(prophageID, MGG_K_locus, phage.type) %>%
    # now make sure all the qulity prophages from this locus are there
    right_join(phages.this.k.locus, by = c("prophageID", "MGG_K_locus"))
  quality.phage.types = rbind(quality.phage.types, phage.types.this.k.locus.at.rbp)
}
quality.phage.types$MGG_K_locus = factor(quality.phage.types$MGG_K_locus, levels = k.locus.order)
quality.phage.types$phage.type = factor(quality.phage.types$phage.type, 
                                            levels = c(top.ecods, "Multiple", "None of above", "No RBP detected"))


# Plot the data
plot.phage.types.at.rbp = ggplot(quality.phage.types %>%
           group_by(MGG_K_locus, phage.type) %>%
           summarise(Num = n_distinct(prophageID)) %>%
           ungroup() %>%
           mutate(phage.type = ECOD_labeller(phage.type))) + 
  geom_col(aes(x = MGG_K_locus, y=Num, fill = phage.type), width = 0.7) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = TEXT_SIZE),
        axis.text.y = element_text(size = TEXT_SIZE),
        axis.title.y = element_text(size = TEXT_SIZE),
        legend.text = element_text(size = TEXT_SIZE),
        legend.title = element_text(size = TEXT_SIZE)) +
  scale_fill_manual(values = custom_colors, name = "ECOD Domain at RBP") +
  xlab("") +
  ylab("Num. high quality\nnearly complete phages")

fig5 = plot.gen.spec / plot.phage.types.at.rbp
fig5 = fig5 +  plot_annotation(tag_levels = 'A') +  theme(plot.title = element_text(size = 20))
ggsave(filename = 'Figure5.jpg', fig5, width = 12, height = 13)

########################################################################
# Look into quality prophages from our K loci and if they have multiple ecods at their RBPs
dir.create('suppl-figs', showWarnings = F)
top.ecods.in.quality.phage.rbps = prophage.proteins.with.ecods %>% 
  inner_join(quality.phage.types) %>%
  filter(PC %in% all.rbp.associated.pcs.unique)  %>% 
  distinct(prophageID, ecod.t) %>%
  filter(ecod.t %in% top.ecods) %>%
  mutate(value = 1) %>%
  tidyr::spread(ecod.t, value, fill = 0) 

jpeg(filename = 'suppl-figs/rbp_top_ecod_combinations.jpg', width = 24, height = 10, units = "cm", res = 300)
upset(data = top.ecods.in.quality.phage.rbps,
      keep.order = TRUE,
      nintersects = 10,
      order.by = "freq",
      nsets = 8, number.angles = 30, point.size = 1.5, line.size = 1, 
      mainbar.y.label = "Num. high quality prophages\nwith combination of ECODs at RBP", sets.x.label = "Num. high quality prophages\nwith this ECOD at RBP", 
      text.scale = c(0.8, 0.6, 0.6, 0.6, 0.6, 1),
      main.bar.color = "darkblue",
      sets.bar.color = "darkblue",
      matrix.color = "darkblue")
# last one is the numbers above bars
# one befpre is names of ECODs
# first is the label
dev.off()


# Look at phages from our K loci that do have RBP but no top ECODs within them, what are neighbouting ECODs
phages.with.no.ecods.at.rbps = prophage.proteins.with.ecods  %>%
  inner_join(quality.phage.types) %>% 
  filter(PC %in% all.rbp.associated.pcs.unique) %>%
  group_by(prophageID) %>%
  summarise(any.top.ecod = any(ecod.t %in% top.ecods))%>%
  ungroup() %>%
  filter(!any.top.ecod) 


ecods.in.phages.with.no.ecods.at.rbps = prophage.proteins.with.ecods %>%
  inner_join(phages.with.no.ecods.at.rbps) %>%
  mutate(is.rbp = PC %in% all.rbp.associated.pcs.unique) %>%
  group_by(prophageID) %>%
  mutate(is.around.rbp = is.near.RBP(is.rbp, protID)) %>%
  ungroup() %>%
  filter(ecod.t != "unknown function" & is.around.rbp) %>%
  group_by(ecod.t) %>%
  summarise(num.phages = n_distinct(prophageID)) %>%
  arrange(desc(num.phages))
ecods.in.phages.with.no.ecods.at.rbps$ecod.t = factor(ecods.in.phages.with.no.ecods.at.rbps$ecod.t,
                                                        levels = unique(ecods.in.phages.with.no.ecods.at.rbps$ecod.t))

ecods.around.phages.with.no.ecod.in.rbp = ggplot(ecods.in.phages.with.no.ecods.at.rbps %>% 
                                                   filter(num.phages > 15)) +
  geom_col(aes(x = ecod.t, y= num.phages), fill = "darkblue") +
  coord_flip() +
  theme_bw() + 
  xlab("ECOD detected near RBP") +
  ylab("Num. high quality prophages with no detected domains at RBP")
ggsave(filename = 'suppl-figs/ecods.around.phages.with.no.ecod.in.rbp.jpg', ecods.around.phages.with.no.ecod.in.rbp, width = 8, height = 8)



################################################ Look at acetyltransferases ################################################
# acetyltrasnsf ECODs in high quality prophages
acetyltransf_ecod_fnames =  ecod_annot_table_f[which(grepl("Acetyltransf", ecod_annot_table_f$ecod.f) | grepl("acetyltransf", ecod_annot_table_f$ecod.f)),] %>% 
  distinct(ecod.t, ecod.f) %>% 
  pull(ecod.f)
# acetyltrasnsf PHROGs in high quality prophages
acetyltransf_phrog_names =  phrog_annot_table[which(grepl("Acetyltransf", phrog_annot_table$PHROG.annotation) | grepl("acetyltransf", phrog_annot_table$PHROG.annotation)),] %>% 
  distinct(PHROG.annotation) %>%
  pull(PHROG.annotation)

quality.phage.proteins = quality_prophages_table %>%
  distinct(prophageID, MGG_K_locus) %>%
  filter(MGG_K_locus %in% k.locus.order) %>%
  left_join(pcs2proteins) 

acet.in.quality.phage.proteins = quality.phage.proteins %>%
  left_join(ecod_annot_table_f, by = "PC", relationship = "many-to-many") %>%
  mutate(ecod.acetyltransferase = ecod.f %in% acetyltransf_ecod_fnames) %>%
  distinct(prophageID, MGG_K_locus, PC, ecod.acetyltransferase) %>%
  left_join(phrog_annot_table, by = "PC", relationship = "many-to-many") %>%
  mutate(phrog.acetyltransferase = PHROG.annotation %in% acetyltransf_phrog_names) %>%
  group_by(MGG_K_locus, prophageID) %>%
  summarise(phrog.acetyl.phage = any(phrog.acetyltransferase) & !any(ecod.acetyltransferase),
            ecod.acetyl.phage = !any(phrog.acetyltransferase) & any(ecod.acetyltransferase),
            phrog.ecod.acetyl.phage = any(phrog.acetyltransferase) & any(ecod.acetyltransferase),
            no.acetyl.phage = !any(phrog.acetyltransferase) & !any(ecod.acetyltransferase)) %>%
  group_by(MGG_K_locus) %>%
  summarise(num.phages = n_distinct(prophageID),
            `PHROG & ECOD` = n_distinct(prophageID[phrog.ecod.acetyl.phage])/num.phages,
            `PHROG` = n_distinct(prophageID[phrog.acetyl.phage])/num.phages,
            `ECOD` = n_distinct(prophageID[ecod.acetyl.phage])/num.phages,
            NEITHER = n_distinct(prophageID[no.acetyl.phage])/num.phages) %>%
  ungroup() %>%
  tidyr::gather(key = "DB", value = prop.phages, `PHROG & ECOD`, PHROG, ECOD, NEITHER)


acet.in.quality.phage.proteins$MGG_K_locus = factor(acet.in.quality.phage.proteins$MGG_K_locus, levels = k.locus.order)
acet.in.quality.phage.proteins$DB = factor(acet.in.quality.phage.proteins$DB, levels = c( "NEITHER", "ECOD", "PHROG", "PHROG & ECOD"))

aceltytransferase_phages = ggplot(acet.in.quality.phage.proteins) + 
  geom_col(aes(x = MGG_K_locus, y=prop.phages, fill = DB)) +
  theme_bw() + 
  scale_fill_manual(values = c('PHROG & ECOD' = "red", "PHROG" = "orange", ECOD = "beige", NEITHER = "gray80"), 
                    name = "Acetyltransferase detected\nin prophage by:") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 10),
        legend.text=element_text(size=10)) +
  xlab("") +
  #ylab("# High quality nearly complete phages")
  ylab("% High quality nearly complete phages")
ggsave(filename =  'suppl-figs/acetyltransf_phages.jpg', aceltytransferase_phages, width = 13, height = 6)


