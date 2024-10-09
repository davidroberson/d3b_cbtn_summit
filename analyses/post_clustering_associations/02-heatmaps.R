# script to create plots for each data modality arranged and annotated by clusters
# to determine how the clusters relate to subtypes

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(randomcoloR)
  library(gplots)
  library(ggplotify)
  library(pheatmap)
  library(gridExtra)
})

# parse command line options
option_list <- list(
  make_option(c("--cluster_file"), type = "character", help = "File with multi-modal derived clusters (.tsv)"),
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--rna_file"), type = "character", help = "RNA expression file (.tsv)"),
  make_option(c("--feature_scores_rna"), type = "character", help = "Output for RNA data after running cluster stats (.tsv)"),
  make_option(c("--cnv_file"), type = "character", help = "CNV file (.tsv)"),
  make_option(c("--feature_scores_cnv"), type = "character", help = "Output for CNV data after running cluster stats (.tsv)"),
  make_option(c("--snv_file"), type = "character", help = "SNV file (.tsv)"),
  make_option(c("--feature_scores_snv"), type = "character", help = "Output for SNV data after running cluster stats (.tsv)"),
  make_option(c("--methyl_file"), type = "character", help = "Methylation file (.tsv)"),
  make_option(
    c("--feature_scores_methyl"),
    type = "character",
    help = "Output for Methylation data after running cluster stats (.tsv)"
  ),
  make_option(c("--splice_file"), type = "character", help = "PSI values (.tsv)"),
  make_option(
    c("--feature_scores_splice"),
    type = "character",
    help = "Output for Splice data after running cluster stats (.tsv)"
  ),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

############################# generate feature-level heatmaps ############################

# feature-level heatmaps using top 10 most representative features per modality

# expression
cat('Plotting top intNMF expression features')
feature_scores_rna <- read_tsv(opt$feature_scores_rna)
feature_scores_rna <- feature_scores_rna %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cluster") %>%
  tidyr::gather(key = "feature", value = "value", -c(cluster)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::slice_head(n = 10) %>%
  tidyr::spread(key = "feature", value = "value", fill = 0) %>%
  tibble::column_to_rownames("cluster")
p1 <- pheatmap::pheatmap(
  mat = t(feature_scores_rna),
  fontsize = 5,
  cellwidth = 8,
  cellheight = 5,
  scale = "row",
  angle_col = 0,
  silent = T,
  main = paste0("Expression Data\nTop 10 features per cluster")
)

# methylation
cat('Plotting top intNMF methylation features')
feature_scores_methyl <- read_tsv(opt$feature_scores_methyl)
feature_scores_methyl <- feature_scores_methyl %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cluster") %>%
  tidyr::gather(key = "feature", value = "value", -c(cluster)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::slice_head(n = 10) %>%
  tidyr::spread(key = "feature", value = "value", fill = 0) %>%
  tibble::column_to_rownames("cluster")
p2 <- pheatmap::pheatmap(
  mat = t(feature_scores_methyl),
  fontsize = 5,
  cellwidth = 8,
  cellheight = 5,
  scale = "row",
  angle_col = 0,
  silent = T,
  main = paste0("Methylation Data\nTop 10 features per cluster")
)

# splicing
cat('Plotting top intNMF splicing features')
feature_scores_splice <- read_tsv(opt$feature_scores_splice)
feature_scores_splice <- feature_scores_splice %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cluster") %>%
  tidyr::gather(key = "feature", value = "value", -c(cluster)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(desc(value)) %>%
  dplyr::slice_head(n = 10) %>%
  tidyr::spread(key = "feature", value = "value", fill = 0) %>%
  tibble::column_to_rownames("cluster")
p3 <- pheatmap::pheatmap(
  mat = t(feature_scores_splice),
  fontsize = 5,
  cellwidth = 8,
  cellheight = 5,
  scale = "row",
  angle_col = 0,
  silent = T,
  main = paste0("Splice Data\nTop 10 features per cluster")
)

# combine all in one
pdf(
  file.path(plots_dir, "feature_level_heatmaps.pdf"),
  width = 12,
  height = 23
)
grid.arrange(arrangeGrob(
  grobs = list(p1$gtable, p2$gtable, p3$gtable),
  ncol = 3
))
dev.off()

############################# create annotation for sample-level heatmap #############################

# multi-modal derived clusters
mm_clusters <- read_tsv(file = opt$cluster_file)

# combine Multi-modal clusters with RNA-derived molecular subtypes
histology_file <- opt$histology_file
anno_file_rna <- read_tsv(file = histology_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype, CNS_region) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  dplyr::inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# combine Multi-modal clusters with methylation-derived subclass
anno_file_methyl <- read_tsv(file = opt$histology_file) %>%
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

############################# generate sample-level heatmaps ############################

# heatmaps ordered by sample cluster, by molecular subtype, methylation derived subclass
# define distinct color palette for each annotation
set.seed(100)
cluster_palette <- distinctColorPalette(length(unique(anno_file$mm_cluster)))
names(cluster_palette) <- sort(unique(anno_file$mm_cluster))

set.seed(100)
mol_subtype_palette <- distinctColorPalette(length(unique(anno_file$molecular_subtype)))
names(mol_subtype_palette) <- sort(unique(anno_file$molecular_subtype))

set.seed(100)
dkfz_v11_methylation_subclass_palette <- distinctColorPalette(length(unique(
  anno_file$dkfz_v11_methylation_subclass
)))
names(dkfz_v11_methylation_subclass_palette) <- sort(unique(anno_file$dkfz_v11_methylation_subclass))

set.seed(100)
dkfz_v12_methylation_subclass_palette <- distinctColorPalette(length(unique(
  anno_file$dkfz_v12_methylation_subclass
)))
names(dkfz_v12_methylation_subclass_palette) <- sort(unique(anno_file$dkfz_v12_methylation_subclass))

# list of annotation and their colors
annots_colors <- list(
  cluster = cluster_palette,
  molecular_subtype = mol_subtype_palette,
  dkfz_v11_methylation_subclass = dkfz_v11_methylation_subclass_palette,
  dkfz_v12_methylation_subclass = dkfz_v12_methylation_subclass_palette
)

# arrange
annots <- anno_file %>%
  tibble::column_to_rownames('sample_id') %>%
  dplyr::select(
    dkfz_v12_methylation_subclass,
    dkfz_v11_methylation_subclass,
    molecular_subtype,
    mm_cluster
  ) %>%
  dplyr::arrange(
    mm_cluster,
    molecular_subtype,
    dkfz_v11_methylation_subclass,
    dkfz_v12_methylation_subclass
  )
annots$mm_cluster <- as.character(annots$mm_cluster)

# expression
pdf(
  file.path(plots_dir, "sample_level_heatmaps.pdf"),
  width = 11,
  height = 12,
  onefile = T
)
count_data <- read_tsv(opt$rna_file) %>% column_to_rownames()
p1 <- count_data[rownames(annots), ] %>%
  dplyr::select(colnames(feature_scores_rna)) %>%
  t() %>%
  pheatmap::pheatmap(
    annotation_col = annots,
    annotation_colors = annots_colors,
    fontsize = 5,
    cellwidth = 2,
    cellheight = 5,
    show_colnames = F,
    cluster_cols = F,
    cluster_rows = T,
    scale = "row",
    color = bluered(256),
    silent = T,
    main = paste0("Expression Data\nTop 10 features per cluster")
  )
p1 <- p1 %>% as.ggplot
print(p1)

# methylation
methyl_data <- read_tsv(opt$methyl_file) %>% column_to_rownames()
p2 <- methyl_data[rownames(annots), ] %>%
  dplyr::select(colnames(feature_scores_methyl)) %>%
  t() %>%
  pheatmap::pheatmap(
    annotation_col = annots,
    annotation_colors = annots_colors,
    fontsize = 5,
    cellwidth = 2,
    cellheight = 5,
    show_colnames = F,
    cluster_cols = F,
    cluster_rows = T,
    scale = "row",
    color = bluered(256),
    silent = T,
    main = paste0("Methylation Data\nTop 10 features per cluster")
  )
p2 <- p2 %>% as.ggplot
print(p2)

# splicing
splice_data <- read_tsv(opt$splice_file) %>% column_to_rownames()
p3 <- splice_data[rownames(annots), ] %>%
  dplyr::select(colnames(feature_scores_splice)) %>%
  t() %>%
  pheatmap::pheatmap(
    annotation_col = annots,
    annotation_colors = annots_colors,
    fontsize = 5,
    cellwidth = 2,
    cellheight = 5,
    show_colnames = F,
    cluster_cols = F,
    cluster_rows = T,
    scale = "row",
    color = bluered(256),
    silent = T,
    main = paste0("Splice Data\nTop 10 features per cluster")
  )
p3 <- p3 %>% as.ggplot
print(p3)
dev.off()
