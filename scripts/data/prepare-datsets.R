### load libraries
suppressWarnings(library(yaml))
suppressWarnings(library(data.table))
suppressWarnings(library(ape))


### load config
config.filename <- file.path(rprojroot::find_rstudio_root_file(), "config", "config.yaml")
cfg <- yaml::read_yaml(config.filename)

####################################################################
############################### MAIN ###############################
####################################################################

# set paths
data.input.path <- file.path(cfg$paths$rafal$main, cfg$paths$mgg$db_input_rel)
data.output.path <- file.path(cfg$paths$rafal$main, cfg$paths$rafal$db_output_rel)
dir.create(data.output.path, recursive = T, showWarnings = F)

# set file paths
bacterial.tree.filename <- file.path(data.input.path, "bacteria_iqtree.nwk")
bacteria.metadata.filename <- file.path(data.input.path, "bacteria.tsv")
prophage.metadata.filename <- file.path(data.input.path, 'prophages.tsv')
gwas.clusters.functions.filename <- file.path(data.input.path, 'clusters_functions_best_all.tsv')
gwas.hits.all.filename <- file.path(data.input.path, 'pyseer_hits_all.tsv')
protein.table.filename <- file.path(data.input.path, 'predictions_and_enzymes.tsv')

bacterial.tree.outfile <- file.path(data.output.path, "bacteria_iqtree.nwk")
bacteria.metadata.outfile <- file.path(data.output.path, "bacteria-metadata.csv")
prophage.metadata.outfile <- file.path(data.output.path, 'prophages-metadata.csv')
gwas.clusters.functions.outfile <- file.path(data.output.path, 'clusters_functions_best_all.csv')
gwas.hits.all.outfile <- file.path(data.output.path, 'pyseer_hits_all.csv')
protein.table.outfile <- file.path(data.output.path, 'predictions_and_enzymes.csv')

# read data
bacterial.tree <- read.tree(bacterial.tree.filename)
bacteria.metadata <- fread(bacteria.metadata.filename)
prophage.metadata <- fread(prophage.metadata.filename)
gwas.clusters.functions <- fread(gwas.clusters.functions.filename)
gwas.hits.all <- fread(gwas.hits.all.filename)
protein.table <- fread(protein.table.filename)

# write data
write.tree(bacterial.tree, file = bacterial.tree.outfile)
fwrite(bacteria.metadata, file = bacteria.metadata.outfile)
fwrite(prophage.metadata, file = prophage.metadata.outfile)
fwrite(gwas.clusters.functions, file = gwas.clusters.functions.outfile)
fwrite(gwas.hits.all, file = gwas.hits.all.outfile)
fwrite(protein.table, file = protein.table.outfile)








