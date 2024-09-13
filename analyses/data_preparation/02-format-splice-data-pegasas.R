# this script is to reduce the size of input splice file for PEGASAS analysis
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

# parse command line options
option_list <- list(
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--short_histology"), type = "character", help = "Short histology of interest"),
  make_option(c("--rmats_splice_file"), type = "character", help = "RMATs file (tsv.gz)"),
  make_option(
    c("--functional_sites_splice_file"),
    type = "character",
    help = "file containing functional sites"
  ),
  make_option(c("--output_dir"), type = "character", help = "path to results directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# subset histology to histology of interest
histology_file <- read_tsv(opt$histology_file)
histology_file <- histology_file %>%
  filter(short_histology == opt$short_histology,
         experimental_strategy == "RNA-Seq")

# read rmats file
splice_mat <- data.table::fread(opt$rmats_splice_file)

# subset splice file to biospecimens of interest
splice_mat <- splice_mat %>%
  filter(sample_id %in% histology_file$Kids_First_Biospecimen_ID)

# subset to SE events only
splice_mat <- splice_mat %>%
  filter(splicing_case == "SE")

# 14005221 rows
print(nrow(splice_mat))

# format per PEGASAS specifications
splice_mat <- splice_mat %>%
  mutate(
    exonStart = exonStart_0base + 1,
    splice_id = paste0(
      geneSymbol,
      "_",
      exonStart,
      "-",
      exonEnd,
      "_",
      upstreamES,
      "-",
      upstreamEE,
      "_",
      downstreamES,
      "-",
      downstreamEE
    )
  )

# filter to functional sites (from Ammar Naqvi using v12)
functional_sites_splice_file <- readRDS(opt$functional_sites_splice_file)
splice_mat <- splice_mat %>%
  filter(splice_id %in% functional_sites_splice_file$Splice_ID)
splice_mat <- splice_mat %>%
  dplyr::select(-c(splice_id))

# 4702509 rows
print(nrow(splice_mat))

# write output file
data.table::fwrite(
  x = splice_mat,
  file = file.path(output_dir, "splice-events-rmats-functional-sites.tsv.gz"),
  sep = "\t",
  na = NA,
  quote = FALSE
)
