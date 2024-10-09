# script to identify differentially methylation CpG sites using limma and pathway enrichment using missMethyl::gsameth

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(msigdbr)
  library(missMethyl)
})

# parse command line options
option_list <- list(
  make_option(c("--methyl_mat"), type = "character", help = "methylation data matrix, preferably m-values (.rds) "),
  make_option(c("--diffexpr_sites"), type = "character", help = "methylation data matrix, preferably m-values (.rds) "),
  make_option(c("--msigdb"), type = "character", help = "reactome, kegg, or hallmark"),
  make_option(c("--prefix"), type = "character", help = "prefix for output files"),
  make_option(c("--output_dir"), type = "character", help = "path to results directory"),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
diffexpr_sites <- opt$diffexpr_sites
methyl_mat <- opt$methyl_mat
msigdb <- opt$msigdb
prefix <- opt$prefix

# output directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# genesets
cat('Reading pathway annotations')
if (msigdb == "kegg") {
  category <- "C2"
  subcategory = "CP:KEGG"
} else if (msigdb == "reactome") {
  category <- "C2"
  subcategory = "CP:REACTOME"
} else if (msigdb == "hallmark") {
  category <- "H"
  subcategory <- NULL
}
gene_set <- msigdbr::msigdbr(species = "Homo sapiens",
                             category = category,
                             subcategory = subcategory)
gene_set <- gene_set %>% dplyr::select(gs_name, entrez_gene)
gene_set_entrez <- base::split(gene_set$entrez_gene, list(gene_set$gs_name))

# read differentially methylated sites
diffexpr_sites <- read_tsv(diffexpr_sites)

# read methylation matrix for pulling all probes used in the analysis
cat('Reading methylation values')
methyl_mat <- readRDS(methyl_mat) %>% 
  dplyr::slice_head(n = 500000) %>%
  na.omit() %>%
  dplyr::filter(!duplicated(Probe_ID))
methyl_mat <- methyl_mat %>%
  tibble::column_to_rownames('Probe_ID')
all_cpgs <- rownames(methyl_mat)

# use a for-loop
clusters <- unique(diffexpr_sites$cluster)
all_pathway_df <- data.frame()
pdf(
  file = file.path(plots_dir, paste0(prefix, "_gsameth_pathways.pdf")),
  width = 12,
  height = 10
)
for (i in 1:length(clusters)) {
  sigCpGs <- diffexpr_sites %>%
    dplyr::filter(cluster == clusters[i])
  
  # pathway enrichment using gsameth
  # take top 10000 instead of all significant probes because the latter results in too many probes and hence pathway enrichment does not work
  sigCpGs <- sigCpGs %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = 10000) %>%
    dplyr::pull(probes)
  pathway_df <- gsameth(
    sig.cpg = sigCpGs,
    all.cpg = all_cpgs,
    collection = gene_set_entrez,
    array.type = "EPIC",
    sig.genes = TRUE
  )
  
  # significant pathways
  pathway_df <- pathway_df %>%
    as.data.frame() %>%
    dplyr::mutate(cluster = clusters[i]) %>%
    tibble::rownames_to_column("pathway") %>%
    dplyr::arrange(FDR) %>%
    dplyr::filter(FDR < 0.1) # not many significant pathways, so use a loose cutoff
  
  # combine with full output
  all_pathway_df <- rbind(all_pathway_df, pathway_df)
  
  # barplot per cluster (plot only top 50)
  if (nrow(pathway_df) > 0) {
    pathway_df <- pathway_df %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 50)
    pathway_df$pathway <- gsub("_", " ", pathway_df$pathway)
    pathway_df$pathway <- factor(pathway_df$pathway, levels = pathway_df$pathway)
    p <- ggplot(pathway_df, aes(pathway, y = (-1) * log10(FDR), fill = FDR)) +
      geom_bar(stat = "identity") + coord_flip() + theme_bw() +
      xlab("") +
      ylab("-log10 (FDR)") +
      theme(plot.margin = unit(c(1, 1, 1, 1), "cm")) +
      scale_x_discrete(
        labels = function(x)
          str_wrap(x, width = 80)
      ) +
      ggtitle(paste(
        "Enrichment for Cluster:",
        i,
        "\nTop 50 Differential Pathways"
      ))
    print(p)
  }
}
dev.off()

# write pathway output to tsv
write_tsv(x = all_pathway_df, file = file.path(
  output_dir,
  paste0(prefix, "_gsameth_output_per_cluster.tsv")
))
