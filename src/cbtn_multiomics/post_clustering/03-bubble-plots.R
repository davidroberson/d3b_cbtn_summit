# script to create plots for each data modality arranged and annotated by clusters
# to determine how the clusters relate to subtypes

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(gplots)
  library(corrplot)
})

# parse command line options
option_list <- list(
  make_option(c("--cluster_file"), type = "character", help = "File with multi-modal derived clusters (.tsv)"),
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# multi-modal derived clusters
mm_clusters <- read_tsv(file = opt$cluster_file)

# combine Multi-modal clusters with RNA-derived molecular subtypes
cat('Joining multi-omic cluster annotations with histologic data \n')
histology_file <- opt$histology_file
anno_file_rna <- read_tsv(file = histology_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype, CNS_region) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  dplyr::inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# combine Multi-modal clusters with methylation-derived subclass
anno_file_methyl <- read_tsv(file = histology_file) %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    dkfz_v11_methylation_subclass,
    dkfz_v12_methylation_subclass
  ) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  dplyr::inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_Methyl")

# combine both and create one standardized annotation file
anno_file <- anno_file_rna %>%
  dplyr::inner_join(anno_file_methyl)
anno_file$dkfz_v11_methylation_subclass <- gsub(", | ", "_", anno_file$dkfz_v11_methylation_subclass)

############################# generate balloon and corrplots ############################

# generate balloon plot with at least 5 samples in a group
# Multi-modal clusters vs RNA-derived molecular subtypes
cat('Plotting multi-omic clusters in relation to known subtypes \n')
dat <- anno_file %>%
  filter(!is.na(molecular_subtype)) %>%
  group_by(molecular_subtype, mm_cluster)  %>%
  summarise(n = n()) %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = mm_cluster, value = n, fill = 0) %>%
  column_to_rownames("molecular_subtype")
pdf(
  file = file.path(plots_dir, "mm_clusters_vs_molsubtype_balloonplot.pdf"),
  width = 10
)
balloonplot(
  x = as.table(as.matrix(t(dat))),
  main = "Multi-modal clusters vs Molecular subtypes",
  xlab = "",
  ylab = "",
  label = TRUE,
  show.margins = FALSE
)
dev.off()

# 4) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# Multi-modal clusters vs RNA-derived molecular subtypes
chisq <- chisq.test(dat)
pdf(file = file.path(plots_dir, "mm_clusters_vs_molsubtype_corrplot.pdf"))
corrplot(
  chisq$residuals,
  is.cor = FALSE,
  tl.srt = 360,
  tl.offset = 1,
  mar = c(1, 2, 1, 1),
  title = "Multi-modal clusters vs Molecular subtypes"
)
dev.off()

# 3) generate balloon plot with at least 5 samples in a group
# Multi-modal clusters vs methylation-derived dkfz_v11_methylation_subclass
dat <- anno_file %>%
  dplyr::filter(!is.na(dkfz_v11_methylation_subclass)) %>%
  dplyr::group_by(dkfz_v11_methylation_subclass, mm_cluster)  %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(nmax = max(n)) %>%
  dplyr::filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  tidyr::spread(key = mm_cluster, value = n, fill = 0) %>%
  tibble::column_to_rownames("dkfz_v11_methylation_subclass")
pdf(
  file = file.path(
    plots_dir,
    "mm_clusters_vs_dkfz_v11_methylation_subclass_balloonplot.pdf"
  ),
  width = 14
)
balloonplot(
  x = as.table(as.matrix(t(dat))),
  main = "Multi-modal clusters vs dkfz_v11_methylation_subclass",
  xlab = "",
  ylab = "",
  label = TRUE,
  show.margins = FALSE
)
dev.off()

# 4) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# Multi-modal clusters vs methylation-derived dkfz_v11_methylation_subclass
chisq <- chisq.test(dat)
pdf(file = file.path(
  plots_dir,
  "mm_clusters_vs_dkfz_v11_methylation_subclass_corrplot.pdf"
))
corrplot(
  chisq$residuals,
  is.cor = FALSE,
  tl.srt = 360,
  tl.offset = 1,
  mar = c(1, 2, 1, 1),
  title = "Multi-modal clusters vs dkfz_v11_methylation_subclass"
)
dev.off()

# 3) generate balloon plot with at least 5 samples in a group
# Multi-modal clusters vs methylation-derived dkfz_v12_methylation_subclass
dat <- anno_file %>%
  dplyr::filter(!is.na(dkfz_v12_methylation_subclass)) %>%
  dplyr::group_by(dkfz_v12_methylation_subclass, mm_cluster)  %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(nmax = max(n)) %>%
  dplyr::filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  tidyr::spread(key = mm_cluster, value = n, fill = 0) %>%
  tibble::column_to_rownames("dkfz_v12_methylation_subclass")
pdf(
  file = file.path(
    plots_dir,
    "mm_clusters_vs_dkfz_v12_methylation_subclass_balloonplot.pdf"
  ),
  width = 14
)
balloonplot(
  x = as.table(as.matrix(t(dat))),
  main = "Multi-modal clusters vs dkfz_v12_methylation_subclass",
  xlab = "",
  ylab = "",
  label = TRUE,
  show.margins = FALSE
)
dev.off()

# 4) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# Multi-modal clusters vs methylation-derived dkfz_v12_methylation_subclass
chisq <- chisq.test(dat)
pdf(file = file.path(
  plots_dir,
  "mm_clusters_vs_dkfz_v12_methylation_subclass_corrplot.pdf"
))
corrplot(
  chisq$residuals,
  is.cor = FALSE,
  tl.srt = 360,
  tl.offset = 1,
  mar = c(1, 2, 1, 1),
  title = "Multi-modal clusters vs dkfz_v12_methylation_subclass"
)
dev.off()
