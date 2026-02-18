### load libraries
suppressWarnings(library(RColorBrewer))
suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))
suppressWarnings(library(ggplot2))
suppressWarnings(library(readxl))


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

# LOAD DATA

data.fn <- file.path(cfg$paths$rafal$main, cfg$paths$rafal$data_s3)
data.abs <- as.data.table(read_xlsx(data.fn, sheet = 1))
data.pr <- as.data.table(read_xlsx(data.fn, sheet = 2))

#############################################
# First, generate a single category "protein"
#############################################

# time columns (e.g. "0 min", ..., "20 min")
time.cols <- grep(" min$", names(data.abs), value = TRUE)
abs.cols <- grep("absorbance", names(data.pr), value = TRUE)

# define groups of conditions
protein.conditions <- c("KL2", "KL55", "KL111")
ssp.conditions <- c("substrate+KL2", "substrate+KL55", "substrate+KL111") # substrate plus protein (SSP)

# generate a mean value for each biological replicate
data.abs.prot <- data.abs[condition %chin% protein.conditions]
data.abs.prot <- data.abs.prot[, c(
  list(`technical replicate` = NA_integer_),  # optional; we'll set later
  lapply(.SD, function(x) mean(x, na.rm = TRUE))
),
by = .(condition, `biological replicate`),
.SDcols = time.cols]
data.abs.prot[, condition := "protein"]
data.abs.prot[, `biological replicate` := seq_len(.N)]   # 1..6
data.abs.prot[, `technical replicate` := 1L]             # single "collapsed" tech replicate

data.pr.prot <- data.pr[condition %chin% protein.conditions]
data.pr.prot <- data.pr.prot[, c(
  list(`technical replicate` = NA_integer_),  # optional; we'll set later
  lapply(.SD, function(x) mean(x, na.rm = TRUE))
),
by = .(condition, `biological replicate`),
.SDcols = abs.cols]
data.pr.prot[, condition := "protein"]
data.pr.prot[, `biological replicate` := seq_len(.N)]   # 1..6
data.pr.prot[, `technical replicate` := 1L]             # single "collapsed" tech replicate

# merge with the rest
data.abs.others <- data.abs[!condition %chin% protein.conditions]
data.abs.processed <- rbind(data.abs.others, data.abs.prot)
data.abs.processed[, `technical replicate` := NULL]
setnames(data.abs.processed, old = "biological replicate", new = "replicate")

data.pr.others <- data.pr[!condition %chin% protein.conditions]
data.pr.processed <- rbind(data.pr.others, data.pr.prot)
data.pr.processed[, `technical replicate` := NULL]
setnames(data.pr.processed, old = "biological replicate", new = "replicate")

#############################################
# Second, generate a single category "protein"
#############################################

data.abs.processed.substrate <- data.abs.processed[condition == "substrate"]
data.abs.processed.substrate.mean <- data.abs.processed.substrate[, lapply(.SD, mean, na.rm = TRUE),  .SDcols = time.cols]
data.abs.processed <- rbind(data.abs.processed, data.abs.processed.substrate.mean, fill = T)
data.abs.processed[is.na(condition)]$condition <- "substrate_mean"

#############################
# Third, subtract background
#############################

bg <- unlist(data.abs.processed[condition == "substrate_mean", ..time.cols], use.names = TRUE)
data.abs.processed.ssp <- data.abs.processed[condition %chin% ssp.conditions]
data.abs.processed.nonssp <- data.abs.processed[!condition %chin% ssp.conditions]
data.abs.processed.ssp[, (time.cols) := Map(`-`, .SD, bg), .SDcols = time.cols]
data.abs.processed.ssp[, condition := sprintf("%s-nbg", condition)]

data.abs.processed.final <- rbind(data.abs.processed.ssp, data.abs.processed.nonssp)
data.abs.processed.final <- data.abs.processed.final[!condition %chin% c("substrate_mean")]

#############################
# Fourth, prepare for plotting
#############################

### Plot panel A

# Melt the data into a long format
longA <- melt(
  data.abs.processed.final,
  id.vars      = c("condition", "replicate"),
  measure.vars = time.cols,
  variable.name = "time",
  value.name    = "value"
)

# make time numeric (minutes)
longA[, time := as.integer(sub(" min$", "", time))]

# summary table (exactly 4 columns)
summary.table <- longA[
  ,
  .(value_mean = mean(value, na.rm = TRUE),
    value_sd   = sd(value,   na.rm = TRUE)),
  by = .(condition, time)
][order(condition, time)]

#############################
# Five, plot absorbance
#############################

condition_order <- c(
  "substrate",
  "protein",
  "substrate+KL111-nbg",
  "substrate+KL2-nbg",
  "substrate+KL55-nbg"
)

summary.table$condition <- factor(
  summary.table$condition,
  levels = condition_order
)

condition_labels <- c(
  "protein"             = "Protein only",
  "substrate"           = "Substrate only",
  "substrate+KL2-nbg"   = "Substrate + 164_08-KL2*",
  "substrate+KL55-nbg"  = "Substrate + 174_38-KL55*",
  "substrate+KL111-nbg" = "Substrate + 184_43-KL111*"
)

condition_colors <- c(
  "protein"             = rgb(44,46,113, maxColorValue = 255),
  "substrate"           = rgb(191,104,163, maxColorValue = 255),
  "substrate+KL2-nbg"   = rgb(170,101,48, maxColorValue = 255),
  "substrate+KL55-nbg"  = rgb(177,159,77, maxColorValue = 255),
  "substrate+KL111-nbg" = rgb(73,127,70, maxColorValue = 255)
)



brewer.palette <- "Set3"
panelA <- ggplot(summary.table, aes(x = time, y = value_mean, color = condition, fill = condition)) +
  geom_ribbon(aes(ymin = value_mean - value_sd, ymax = value_mean + value_sd), alpha = 0.4, color = "grey", linewidth = .1) +
  geom_line(size = 1) +
  labs(x = "Time (min)",
       y = "Absorbance at 405 nm (mean ± SD)",
       color = "Condition",
       fill = "Condition") +
  scale_color_manual(values = condition_colors,
                     labels = condition_labels) +
  scale_fill_manual(values = condition_colors,
                    labels = condition_labels) +
  scale_y_continuous(limits = c(-0.01, max(summary.table$value_mean)*1.1)) +
  theme_classic(base_size = 13)
panelA

### Plot panel B


# Melt the data into a long format
longB <- melt(
  data.pr.processed,
  id.vars      = c("condition", "replicate"),
  measure.vars = abs.cols,
  value.name    = "value"
)

# summary table (exactly 4 columns)
summary.table.pr <- longB[
  ,
  .(value_mean = mean(value, na.rm = TRUE),
    value_sd   = sd(value,   na.rm = TRUE)),
  by = .(condition)
][order(condition)]

condition_order <- c(
  "substrate",
  "protein",
  "substrate+KL111",
  "substrate+KL2",
  "substrate+KL55"
)

summary.table.pr$condition <- factor(
  summary.table.pr$condition,
  levels = condition_order
)

condition_labels <- c(
  "protein"             = "Protein only",
  "substrate"           = "Substrate only",
  "substrate+KL2"   = "Substrate +\n 164_08-KL2",
  "substrate+KL55"  = "Substrate +\n 174_38-KL55",
  "substrate+KL111" = "Substrate +\n 184_43-KL111"
)

names(condition_colors) <- gsub("-nbg", "", names(condition_colors))

panelB <- ggplot(summary.table.pr, aes(x = condition, y = value_mean, color = condition)) +
  geom_point(size = 3, alpha = 0.99) +
  geom_errorbar(aes(ymin = value_mean - value_sd,
                    ymax = value_mean + value_sd),
                width = 0.2,
                linewidth = 0.6) +
  scale_color_manual(values = condition_colors,
                     labels = condition_labels) +
  scale_x_discrete(labels = condition_labels) +
  labs(x = "",
       y = "Absorbance at 560 nm (mean ± SD)",
       color = "Condition") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Save files

panelA.outfile <- "Fig7A.pdf"
panelB.outfile <- "Fig7B.pdf"
ggsave(panelA.outfile, panelA, width = 6.5, height = 3.6)
ggsave(panelB.outfile, panelB, width = 4.5, height = 4.5)


