# run multi-modal clustering
suppressPackageStartupMessages({
  library(optparse)
  library(IntNMF)
  library(tidyverse)
})

# parse command line options
option_list <- list(
  make_option(c("--count_file"), type = "character", help = "Counts file (.tsv)"),
  make_option(c("--cnv_file"), type = "character", help = "CNV file (.tsv)"),
  make_option(c("--snv_file"), type = "character", help = "SNV file (.tsv)"),
  make_option(c("--methyl_file"), type = "character", help = "Methylation file (.tsv)"),
  make_option(c("--splice_file"), type = "character", help = "PSI values (.tsv)"),
  make_option(c("--samples_map"), type = "character", help = "Mapping file with bs ids and samples ids (.tsv)"),
  make_option(c("--output_dir"), type = "character", help = "path to results directory"),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# output directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# source function
source(file.path("utils", "run_clusterstats.R"))

# read data matrices and add to a list
dat <- list()

count_data <- read_tsv(opt$count_file) %>%
  column_to_rownames() %>%
  as.matrix()
dat[["rna_data"]] <- count_data
methyl_data <- read_tsv(opt$methyl_file) %>%
  column_to_rownames() %>%
  as.matrix()
dat[["methyl_data"]] <- methyl_data
snv_data <- read_tsv(opt$snv_file) %>%
  column_to_rownames() %>%
  as.matrix()
dat[["snv_data"]] <- snv_data
cnv_data <- read_tsv(opt$cnv_file) %>%
  column_to_rownames() %>%
  as.matrix()
dat[["cnv_data"]] <- cnv_data
splice_data <- read_tsv(opt$splice_file) %>%
  column_to_rownames() %>%
  as.matrix()
dat[["splice_data"]] <- splice_data

# add weights to each data modality
wt = if (is.list(dat))
  rep(1, length(dat)) else 1

# run run_clusterstats across all k-values
# get the nmf output corresponding to the most optimal k
nmf_output <- run_clusterstats(
  dat = dat,
  output_dir = output_dir,
  k_value = 15
)

# save feature scores per cluster for downstream processing
subdir <- file.path(output_dir, "feature_scores")
dir.create(subdir, showWarnings = F, recursive = T)
write_tsv(as.data.frame(nmf_output$H$rna_data),
          file = file.path(subdir, "feature_scores_rna.tsv"))
write_tsv(as.data.frame(nmf_output$H$cnv_data),
          file = file.path(subdir, "feature_scores_cnv.tsv"))
write_tsv(
  as.data.frame(nmf_output$H$methyl_data),
  file = file.path(subdir, "feature_scores_methyl.tsv")
)
write_tsv(as.data.frame(nmf_output$H$snv_data),
          file = file.path(subdir, "feature_scores_snv.tsv"))
write_tsv(
  as.data.frame(nmf_output$H$splice_data),
  file = file.path(subdir, "feature_scores_splice.tsv")
)

# 1) ConsensusMatPlot
# Given the integrative NMF fit object, the function creates image plot of the consensus matrix ordered
# according to clusters groups. Cleaner block structure indicates stronger clusters
pdf(
  file = file.path(plots_dir, "intnmf_consensus_plot.pdf"),
  width = 10,
  height = 10,
  onefile = F
)
ConsensusMatPlot(fit = nmf_output,
                 rowLab = TRUE,
                 colLab = TRUE)
dev.off()

# 2) SilhouettePlot
# Silhouette width plot is returned together with mean silhouette width for each group, overall silhouette width and summary statistics.
pdf(
  file = file.path(plots_dir, "intnmf_silhouette_plot.pdf"),
  width = 10,
  height = 10,
  onefile = F
)
SilhouettePlot(fit = nmf_output, cluster.col = NULL)
dev.off()

# 3) output clusters with per sample
df <- data.frame(sample_id = names(nmf_output$clusters),
                 mm_cluster = nmf_output$clusters)

# combine with sample mapping file that has biospecimens for each data modality
samples_map <- read_tsv(opt$samples_map)
df <- samples_map %>%
  inner_join(df)

# write output
write_tsv(df, file = file.path(output_dir, "intnmf_clusters.tsv"))
