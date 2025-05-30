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

ecod.domains.path <- file.path(cfg$paths$rafal$main, cfg$paths$mgg$ecod_domains_input)
ecod.domains <- fread(
  ecod.domains.path,
  skip = 4,
  colClasses = c("character", rep("auto", 15))  # make the first column a character
)
setnames(ecod.domains, old = "#uid", new = "uid")

# Separate 'f_id' into 'x_level', 'h_level', 't_level', 'f_level'
ecod.domains[, c("x_level", "h_level", "t_level", "f_level") := tstrsplit(f_id, ".", fixed = TRUE)]

# Create hierarchical IDs
ecod.domains[, x_id := x_level]
ecod.domains[, h_id := paste(x_level, h_level, sep = ".")]
ecod.domains[, t_id := paste(x_level, h_level, t_level, sep = ".")]
ecod.domains[, f_id := paste(x_level, h_level, t_level, f_level, sep = ".")]

# propagate T names upwards if missing
missing.h <- ecod.domains$h_name == "NO_H_NAME"
missing.x <- ecod.domains$x_name == "NO_X_NAME"
ecod.domains[missing.h]$h_name <- ecod.domains[missing.h]$t_name
ecod.domains[missing.x]$x_name <- ecod.domains[missing.x]$t_name

ecod.domains$full_id <- sprintf("ECOD_%s_%s", ecod.domains$uid, ecod.domains$ecod_domain_id)

ecod.domains.output <- cfg$data$ecod_domains
fwrite(ecod.domains, file = ecod.domains.output)
