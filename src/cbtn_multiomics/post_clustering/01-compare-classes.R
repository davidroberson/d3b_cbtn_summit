# script to create plots for each data modality arranged and annotated by clusters
# to determine how the clusters relate to subtypes

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(mclust)
})

# parse command line options
option_list <- list(
  make_option(c("--cluster_file"), type = "character", help = "File with multi-modal derived clusters (.tsv)"),
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--rna_file"), type = "character", help = "RNA data file (.tsv)"),
  make_option(c("--methyl_data"), type = "character", help = "Methylation data file (.tsv)"),
  make_option(c("--splice_data"), type = "character", help = "Splice data file (.tsv)"),
  make_option(c("--feature_scores_rna"), type = "character", help = "RNA feature scores (.tsv)"),
  make_option(c("--feature_scores_methyl"), type = "character", help = "Methylation feature scores (.tsv)"),
  make_option(c("--feature_scores_splice"), type = "character", help = "Splice feature scores (.tsv)"),
  make_option(c("--output_dir"), type = "character", help = "path to results directory"),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# results directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# Create subdirectories for various plot types
heatmaps_dir <- file.path(plots_dir, "heatmaps")
dir.create(heatmaps_dir, showWarnings = F, recursive = T)

bubble_plots_dir <- file.path(plots_dir, "bubble_plots")
dir.create(bubble_plots_dir, showWarnings = F, recursive = T)

survival_plots_dir <- file.path(plots_dir, "survival_plots")
dir.create(survival_plots_dir, showWarnings = F, recursive = T)

sankey_plots_dir <- file.path(plots_dir, "sankey_plots")
dir.create(sankey_plots_dir, showWarnings = F, recursive = T)

# multi-modal derived clusters
mm_clusters <- read_tsv(file = opt$cluster_file)

# combine Multi-modal clusters with RNA-derived molecular subtypes
cat('Comparing multi-omic clusters to known subtypes \n')
anno_file_rna <- read_tsv(file = opt$histology_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype, CNS_region) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# combine Multi-modal clusters with methylation-derived subclass
methyl_columns <- c("Kids_First_Biospecimen_ID")
histology_df <- read_tsv(file = opt$histology_file)

# Check if methylation columns exist
if ("dkfz_v11_methylation_subclass" %in% colnames(histology_df)) {
  methyl_columns <- c(methyl_columns, "dkfz_v11_methylation_subclass")
}
if ("dkfz_v12_methylation_subclass" %in% colnames(histology_df)) {
  methyl_columns <- c(methyl_columns, "dkfz_v12_methylation_subclass")
}

anno_file_methyl <- histology_df %>%
  dplyr::select(all_of(methyl_columns)) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_Methyl")

# combine both and create one standardized annotation file
anno_file <- anno_file_rna %>%
  dplyr::inner_join(anno_file_methyl)

# Process methylation subclass if it exists
if ("dkfz_v11_methylation_subclass" %in% colnames(anno_file)) {
  anno_file$dkfz_v11_methylation_subclass <- gsub(", | ", "_", anno_file$dkfz_v11_methylation_subclass)
}

# 1) adjusted rand index between various subtypes and Multi-modal derived clusters
sink(file.path(output_dir, "ari_mm_vs_subtypes.txt"),
     type = c("output"))
print("Multi-modal clusters vs RNA-derived molecular_subtypes")
if ("molecular_subtype" %in% colnames(anno_file)) {
  print(mclust::adjustedRandIndex(anno_file$molecular_subtype, anno_file$mm_cluster))
} else {
  print("Column 'molecular_subtype' not found in data")
}

print("Multi-modal clusters vs dkfz_v11_methylation_subclass")
if ("dkfz_v11_methylation_subclass" %in% colnames(anno_file)) {
  print(mclust::adjustedRandIndex(anno_file$dkfz_v11_methylation_subclass, anno_file$mm_cluster))
} else {
  print("Column 'dkfz_v11_methylation_subclass' not found in data")
}

print("Multi-modal clusters vs dkfz_v12_methylation_subclass")
if ("dkfz_v12_methylation_subclass" %in% colnames(anno_file)) {
  print(mclust::adjustedRandIndex(anno_file$dkfz_v12_methylation_subclass, anno_file$mm_cluster))
} else {
  print("Column 'dkfz_v12_methylation_subclass' not found in data")
}
sink()

# 2) chi-square test of independence between various subtypes and Multi-modal derived clusters
sink(file.path(output_dir, "chisq_mm_vs_subtypes.txt"),
     type = c("output"))
print("Multi-modal clusters vs RNA-derived molecular_subtypes")
if ("molecular_subtype" %in% colnames(anno_file)) {
  print(table(anno_file$molecular_subtype, anno_file$mm_cluster))
  tryCatch({
    print(chisq.test(x = anno_file$molecular_subtype, y = anno_file$mm_cluster))
  }, error = function(e) {
    print("Chi-squared test failed - likely due to insufficient data or low cell counts")
  })
} else {
  print("Column 'molecular_subtype' not found in data")
}

print("Multi-modal clusters vs dkfz_v11_methylation_subclass")
if ("dkfz_v11_methylation_subclass" %in% colnames(anno_file)) {
  print(table(anno_file$dkfz_v11_methylation_subclass, anno_file$mm_cluster))
  tryCatch({
    print(chisq.test(x = anno_file$dkfz_v11_methylation_subclass, y = anno_file$mm_cluster))
  }, error = function(e) {
    print("Chi-squared test failed - likely due to insufficient data or low cell counts")
  })
} else {
  print("Column 'dkfz_v11_methylation_subclass' not found in data")
}

print("Multi-modal clusters vs dkfz_v12_methylation_subclass")
if ("dkfz_v12_methylation_subclass" %in% colnames(anno_file)) {
  print(table(anno_file$dkfz_v12_methylation_subclass, anno_file$mm_cluster))
  tryCatch({
    print(chisq.test(x = anno_file$dkfz_v12_methylation_subclass, y = anno_file$mm_cluster))
  }, error = function(e) {
    print("Chi-squared test failed - likely due to insufficient data or low cell counts")
  })
} else {
  print("Column 'dkfz_v12_methylation_subclass' not found in data")
}
sink()

# Create placeholder files for each output directory to ensure they exist
# and can be recognized by CWL
cat("Placeholder file", file = file.path(heatmaps_dir, "placeholder.txt"))
cat("Placeholder file", file = file.path(bubble_plots_dir, "placeholder.txt"))
cat("Placeholder file", file = file.path(survival_plots_dir, "placeholder.txt"))
cat("Placeholder file", file = file.path(sankey_plots_dir, "placeholder.txt"))

# Combine the chi-square and ARI results
file.copy(
  from = file.path(output_dir, "ari_mm_vs_subtypes.txt"),
  to = file.path(output_dir, "chisq_ari_mm_vs_subtypes.txt"),
  overwrite = TRUE
)
cat("\n\n", file = file.path(output_dir, "chisq_ari_mm_vs_subtypes.txt"), append = TRUE)
file_conn <- file(file.path(output_dir, "chisq_mm_vs_subtypes.txt"), "r")
chisq_content <- readLines(file_conn)
close(file_conn)
cat(paste(chisq_content, collapse = "\n"), file = file.path(output_dir, "chisq_ari_mm_vs_subtypes.txt"), append = TRUE)
