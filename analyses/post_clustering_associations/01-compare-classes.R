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
  make_option(c("--output_dir"), type = "character", help = "path to results directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# results directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# multi-modal derived clusters
mm_clusters <- read_tsv(file = opt$cluster_file)

# combine Multi-modal clusters with RNA-derived molecular subtypes
anno_file_rna <- read_tsv(file = opt$histology_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype, CNS_region) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# combine Multi-modal clusters with methylation-derived subclass
anno_file_methyl <- read_tsv(file = opt$histology_file) %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    dkfz_v11_methylation_subclass,
    dkfz_v12_methylation_subclass
  ) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_Methyl")

# combine both and create one standardized annotation file
anno_file <- anno_file_rna %>%
  inner_join(anno_file_methyl)
anno_file$dkfz_v11_methylation_subclass <- gsub(", | ", "_", anno_file$dkfz_v11_methylation_subclass)

# 1) adjusted rand index between various subtypes and Multi-modal derived clusters
sink(file.path(output_dir, "ari_mm_vs_subtypes.txt"),
     type = c("output"))
print("Multi-modal clusters vs RNA-derived molecular_subtypes")
mclust::adjustedRandIndex(anno_file$molecular_subtype, anno_file$mm_cluster)

print("Multi-modal clusters vs dkfz_v11_methylation_subclass")
mclust::adjustedRandIndex(anno_file$dkfz_v11_methylation_subclass,
                          anno_file$mm_cluster)

print("Multi-modal clusters vs dkfz_v12_methylation_subclass")
mclust::adjustedRandIndex(anno_file$dkfz_v12_methylation_subclass,
                          anno_file$mm_cluster)
sink()

# 2) chi-square test of independence between various subtypes and Multi-modal derived clusters
sink(file.path(output_dir, "chisq_mm_vs_subtypes.txt"),
     type = c("output"))
print("Multi-modal clusters vs RNA-derived molecular_subtypes")
table(anno_file$molecular_subtype, anno_file$mm_cluster)
chisq.test(x = anno_file$molecular_subtype, y = anno_file$mm_cluster)

print("Multi-modal clusters vs dkfz_v11_methylation_subclass")
table(anno_file$dkfz_v11_methylation_subclass,
      anno_file$mm_cluster)
chisq.test(x = anno_file$dkfz_v11_methylation_subclass,
           y = anno_file$mm_cluster)

print("Multi-modal clusters vs dkfz_v12_methylation_subclass")
table(anno_file$dkfz_v12_methylation_subclass,
      anno_file$mm_cluster)
chisq.test(x = anno_file$dkfz_v12_methylation_subclass,
           y = anno_file$mm_cluster)
sink()
