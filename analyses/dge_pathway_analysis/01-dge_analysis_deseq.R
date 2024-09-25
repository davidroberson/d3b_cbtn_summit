# Function: DGE analysis by DESeq2

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(DESeq2)
  library(rtracklayer)
  library(msigdbr)
})

# parse command line options
option_list <- list(
  make_option(c("--expr_mat"), type = "character", help = "expression data matrix, preferably counts (.rds) "),
  make_option(c("--gtf_file"), type = "character", help = "gencode gtf file"),
  make_option(c("--cluster_file"), type = "character", help = "path to cluster annotation file"),
  make_option(c("--results_dir"), type = "character", help = "path to results directory"),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# source function for pathway enrichment
source(file.path("utils", "perform_enrichment_gsea.R"))

# results directory
results_dir <- opt$results_dir
dir.create(results_dir, showWarnings = F, recursive = T)
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# read gtf and filter to protein coding
gencode_gtf <- rtracklayer::import(con = opt$gtf_file) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  dplyr::filter(gene_type == "protein_coding") %>%
  unique()

# count data
expr_mat <- readRDS(opt$expr_mat)

# read cluster information
mm_clusters <- read_tsv(file.path(opt$cluster_file))
mm_clusters <- mm_clusters %>%
  dplyr::arrange(mm_cluster) %>%
  dplyr::mutate(mm_cluster = paste0("cluster_", mm_cluster))

# get sample map information to assign sample ids
expr_mat <- expr_mat %>%
  dplyr::select(mm_clusters$Kids_First_Biospecimen_ID_RNA)
colnames(expr_mat) <- mm_clusters$sample_id

# filter expression count file to contain only protein coding gene
expr_mat <- expr_mat %>%
  filter(rownames(expr_mat) %in% gencode_gtf$gene_name)

# match samples in annotation file and expression matrix
expr_mat <- expr_mat %>%
  dplyr::select(mm_clusters$sample_id)
stopifnot(identical(colnames(expr_mat), mm_clusters$sample_id))

# use a for-loop
clusters <- unique(mm_clusters$mm_cluster)
output_df <- data.frame()
for (i in 1:length(clusters)) {
  # set seed for reproducibility
  set.seed(100)
  
  # prefix for output files and plots
  prefix <- paste0(clusters[i], "_vs_rest")
  print(prefix)
  
  # cluster of interest
  mm_clusters$group <- ifelse(mm_clusters$mm_cluster == clusters[i], "COI", "Others")
  group <- as.factor(mm_clusters$group)
  
  # DESeq2 analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(expr_mat),
    colData = mm_clusters,
    design = ~ group
  )
  dds <- DESeq2::DESeq(dds)
  deseq_results <- DESeq2::results(dds, contrast = c("group", "COI", "Others"))
  deseq_output <- deseq_results %>%
    as.data.frame() %>%
    filter(padj < 0.05) %>%
    dplyr::rename("log2FC" = "log2FoldChange") %>%
    mutate(direction = ifelse(log2FC > 0, "up", "down"),
           comparison = prefix) %>%
    rownames_to_column("genes")
  
  # combine with dataframe
  output_df <- rbind(output_df, deseq_output)
  
  # pathway enrichment using REACTOME
  reactome_pathways <- msigdbr::msigdbr(category = "C2", subcategory = "CP:REACTOME")
  reactome_pathways <- reactome_pathways %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename("term" = "gs_name", "gene" = "gene_symbol")
  perform_enrichment_gsea(
    diffexpr_res = deseq_output,
    pathways = reactome_pathways,
    minGSSize = 10,
    maxGSSize = 150,
    prefix = prefix,
    plots_dir = file.path(plots_dir, "reactome"),
    results_dir = file.path(results_dir, "reactome")
  )
  
  # pathway enrichment using HALLMARK
  hallmark_pathways <- msigdbr::msigdbr(category = "H", subcategory = NULL)
  hallmark_pathways <- hallmark_pathways %>%
    dplyr::select(gs_name, gene_symbol) %>%
    dplyr::rename("term" = "gs_name", "gene" = "gene_symbol")
  perform_enrichment_gsea(
    diffexpr_res = deseq_output,
    pathways = hallmark_pathways,
    minGSSize = 10,
    maxGSSize = 150,
    prefix = prefix,
    plots_dir = file.path(plots_dir, "hallmark"),
    results_dir = file.path(results_dir, "hallmark")
  )
}

# write output to tsv
write_tsv(x = output_df,
          file = file.path(results_dir, "diffexpr_output_per_cluster.tsv"))
