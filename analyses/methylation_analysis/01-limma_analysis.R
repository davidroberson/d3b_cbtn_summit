# script to identify differentially methylation CpG sites using limma

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(limma)
  library(optparse)
  library(data.table)
  library(dplyr)
})

# parse command line options
option_list <- list(
  make_option(c("--methyl_mat"), type = "character", help = "methylation data matrix, preferably m-values (.rds) "),
  make_option(c("--methyl_annot"), type = "character", help = "methylation annotation file (.tsv.gz) "),
  make_option(c("--cluster_file"), type = "character", help = "path to cluster annotation file"),
  make_option(c("--output_dir"), type = "character", help = "path to results directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# output directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# read m-values
methyl_m_values_full <- readRDS(opt$methyl_mat) %>% 
  dplyr::slice_head(n = 500000) %>%
  na.omit() %>%
  dplyr::filter(!duplicated(Probe_ID))
methyl_m_values_full <- methyl_m_values_full %>%
  tibble::column_to_rownames('Probe_ID')

# read annotation
methyl_annot_full <- data.table::fread(opt$methyl_annot) %>%
  dplyr::filter(!duplicated(Probe_ID))

# create generalized function for gene feature analysis
run_analysis <- function(methyl_m_values_full,
                         methyl_annot_full,
                         gene_feature_filter,
                         prefix) {
  # filter only to gene feature using probe annotation file
  methyl_annot <- methyl_annot_full %>%
    dplyr::filter(Gene_Feature %in% gene_feature_filter) %>%
    unique()
  methyl_m_values <- methyl_m_values_full %>%
    dplyr::filter(rownames(methyl_m_values_full) %in% methyl_annot$Probe_ID)
  
  # read cluster information for these samples
  mm_clusters <- read_tsv(file.path(opt$cluster_file))
  mm_clusters <- mm_clusters %>%
    dplyr::arrange(mm_cluster) %>%
    dplyr::mutate(mm_cluster = paste0("cluster_", mm_cluster))
  
  # assign sample ids
  methyl_m_values <- methyl_m_values %>%
    dplyr::select(mm_clusters$Kids_First_Biospecimen_ID_Methyl)
  stopifnot(identical(
    colnames(methyl_m_values),
    mm_clusters$Kids_First_Biospecimen_ID_Methyl
  ))
  colnames(methyl_m_values) <- mm_clusters$sample_id
  
  # match matrix to annotation
  methyl_m_values <- methyl_m_values %>%
    dplyr::select(mm_clusters$sample_id)
  stopifnot(identical(colnames(methyl_m_values), mm_clusters$sample_id))
  
  # use a for-loop
  clusters <- unique(mm_clusters$mm_cluster)
  sigCpGs_output_df <- data.frame()
  for (i in 1:length(clusters)) {
    # cluster of interest
    mm_clusters$group <- ifelse(mm_clusters$mm_cluster == clusters[i], "COI", "Others")
    
    # create design
    mm_clusters$group <- factor(mm_clusters$group, levels = c("Others", "COI"))
    group <- mm_clusters$group
    design <- model.matrix( ~ group)
    rownames(design) <- mm_clusters$sample_id
    
    # fit
    fit.reduced <- limma::lmFit(methyl_m_values, design)
    fit.reduced <- limma::eBayes(fit.reduced, robust = TRUE)
    
    # differentially expressed probes
    toptable_output <- limma::topTable(fit.reduced, coef = 2, n = Inf)
    sigCpGs <- toptable_output %>%
      dplyr::filter(adj.P.Val < 0.05) %>%
      dplyr::mutate(cluster = clusters[i]) %>%
      rownames_to_column("probes")
    print(head(sigCpGs))
    
    # combine with full output
    sigCpGs_output_df <- rbind(sigCpGs_output_df, sigCpGs)
  }
  
  # write significant probes to tsv
  write_tsv(x = sigCpGs_output_df, file = file.path(
    output_dir,
    paste0(prefix, "_diffexpr_probes_per_cluster.tsv")
  ))
}

# run promoter analysis
run_analysis(
  methyl_m_values_full,
  methyl_annot_full,
  gene_feature_filter = "promoter",
  prefix = "promoter"
)

# run gene body analysis
run_analysis(
  methyl_m_values_full,
  methyl_annot_full,
  gene_feature_filter = c("exon", "intron"),
  prefix = "gene_body"
)

# run both gene body + promoter analysis
run_analysis(
  methyl_m_values_full,
  methyl_annot_full,
  gene_feature_filter = c("exon", "intron", "promoter"),
  prefix = "genebody_promoter"
)
