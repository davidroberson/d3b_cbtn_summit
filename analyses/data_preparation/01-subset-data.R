# script to create subset matrices containing only histology of interest

# load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

# parse command line options
option_list <- list(
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--short_histology"), type = "character", help = "Short histology of interest"),
  make_option(c("--count_file"), type = "character", help = "counts matrix (.rds)"),
  make_option(c("--tpm_file"), type = "character", help = "TPM matrix (.rds)"),
  make_option(c("--snv_file"), type = "character", help = "snv consensus + hotspots (.tsv.gz)"),
  make_option(c("--methyl_m_file"), type = "character", help = "methylation m-values (.rds)"),
  make_option(c("--methyl_b_file"), type = "character", help = "methylation beta-values (.rds)"),
  make_option(c("--splice_file"), type = "character", help = "splice file (.rds)"),
  make_option(c("--cnv_file"), type = "character", help = "cnv gainloss.txt (.txt)"),
  make_option(c("--output_dir"), type = "character", help = "path to results directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# read histology and subset to short histology of interest
histology_file <- read_tsv(opt$histology_file) %>%
  dplyr::filter(short_histology %in% opt$short_histology)
write_tsv(histology_file, file = file.path(output_dir, paste0(opt$short_histology, '_', "histologies.tsv")))

# RNA-seq counts
rna_seq_counts <- readRDS(opt$count_file)
rna_seq_counts <- rna_seq_counts %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
saveRDS(
  rna_seq_counts,
  file = file.path(output_dir, paste0(opt$short_histology, '_', "gene-counts-rsem-expected_count-collapsed.rds"))
)
rm(rna_seq_counts)

# RNA-seq TPM
rna_seq_tpm <- readRDS(opt$tpm_file)
rna_seq_tpm <- rna_seq_tpm %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
saveRDS(rna_seq_tpm,
        file = file.path(output_dir, paste0(opt$short_histology, '_', "gene-expression-rsem-tpm-collapsed.rds")))
rm(rna_seq_tpm)

# SNV consensus + hotspots
snv_data <- data.table::fread(opt$snv_file)
snv_data <- snv_data %>%
  dplyr::filter(Tumor_Sample_Barcode %in% histology_file$Kids_First_Biospecimen_ID)
data.table::fwrite(
  snv_data,
  file = file.path(output_dir, paste0(opt$short_histology, '_', "snv-consensus-plus-hotspots.maf.tsv.gz")),
  sep = "\t",
  na = NA,
  quote = FALSE
)
rm(snv_data)

# Methylation m-values
methyl_m_values <- readRDS(opt$methyl_m_file)
methyl_m_values <- methyl_m_values %>%
  dplyr::mutate(Probe_ID = make.names(methyl_m_values$Probe_ID, unique = TRUE)) %>%
  tibble::column_to_rownames("Probe_ID") %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
saveRDS(methyl_m_values, file = file.path(output_dir, paste0(opt$short_histology, '_', "methyl-m-values.rds")))
rm(methyl_m_values)

# Methylation beta-values
methyl_beta_values <- readRDS(opt$methyl_b_file)
methyl_beta_values <- methyl_beta_values %>%
  dplyr::mutate(Probe_ID = make.names(methyl_beta_values$Probe_ID, unique = TRUE)) %>%
  tibble::column_to_rownames("Probe_ID") %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
saveRDS(methyl_beta_values,
        file = file.path(output_dir, paste0(opt$short_histology, '_',"methyl-beta-values.rds")))
rm(methyl_beta_values)

# Splice data (filtered to functional sites from Ammar Naqvi)
splice_data <- readRDS(opt$splice_file)
splice_data <- splice_data %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
saveRDS(splice_data,
        file = file.path(
          output_dir,
          paste0(opt$short_histology, '_', "splice_events_pan_cancer_functional_filter.rds"
        )))
rm(splice_data)

# CNV gainloss file
cnv_data <- data.table::fread(opt$cnv_file)
cnv_data <- cnv_data %>%
  dplyr::filter(BS_ID %in% histology_file$Kids_First_Biospecimen_ID)
data.table::fwrite(
  cnv_data,
  file = file.path(output_dir, paste0(opt$short_histology, '_', "All.gainloss.txt.gz"),
  sep = "\t",
  na = NA,
  quote = FALSE
))
rm(cnv_data)
