# script to identify differentially methylation regions using DMRcate::dmrcate and pathway enrichment using missMethyl::gsaregion

suppressPackageStartupMessages({
  library(tidyverse)
  library(msigdbr)
  library(limma)
  library(DMRcate)
  library(missMethyl)
  library(optparse)
})

mem.maxVSize(vsize = 102400)

# parse command line options
option_list <- list(
  make_option(c("--methyl_mat"), type = "character", help = "methylation data matrix, preferably m-values (.rds) "),
  make_option(c("--methyl_annot"), type = "character", help = "methylation annotation file (.tsv.gz) "),
  make_option(c("--cluster_file"), type = "character", help = "path to cluster annotation file"),
  make_option(c("--msigdb"), type = "character", help = "reactome, kegg, or hallmark"),
  make_option(c("--output_dir"), type = "character", help = "path to results directory"),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
msigdb <- opt$msigdb

# output directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# genesets
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

# read m-values
methyl_m_values_full <- readRDS(opt$methyl_mat) %>% 
  dplyr::slice_head(n = 500000) %>%
  na.omit() %>%
  dplyr::filter(!duplicated(Probe_ID))
methyl_m_values_full <- methyl_m_values_full %>%
  tibble::column_to_rownames("Probe_ID")

# read annotation
methyl_annot_full <- data.table::fread(opt$methyl_annot)

# create generalized function for gene feature analysis
run_analysis <- function(methyl_m_values_full,
                         methyl_annot_full,
                         gene_feature_filter,
                         prefix) {
  # filter only to gene feature using probe annotation file
  methyl_annot <- methyl_annot_full %>%
    filter(Gene_Feature %in% gene_feature_filter) %>%
    unique()
  methyl_m_values <- methyl_m_values_full %>%
    filter(rownames(methyl_m_values_full) %in% methyl_annot$Probe_ID)
  
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
  output_df <- data.frame()
  pdf(
    file = file.path(plots_dir, paste0(prefix, "_gsaregion_pathways.pdf")),
    width = 12,
    height = 10
  )
  for (i in 1:length(clusters)) {
    # cluster of interest
    mm_clusters$group <- ifelse(mm_clusters$mm_cluster == clusters[i], "COI", "Others")
    
    # create design
    mm_clusters$group <- factor(mm_clusters$group, levels = c("Others", "COI"))
    group <- mm_clusters$group
    design <- model.matrix( ~ group)
    rownames(design) <- mm_clusters$sample_id
    
    # run
    myAnnotation <- cpg.annotate(
      object = as.matrix(methyl_m_values),
      datatype = "array",
      what = "M",
      arraytype = "EPICv1",
      analysis.type = "differential",
      design = design,
      coef = 2
    )
    DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2)
    results.ranges <- extractRanges(DMRs)
    
    # gsaregion analysis
    gsa.region <- gsaregion(
      regions = results.ranges,
      all.cpg = rownames(methyl_m_values),
      collection = gene_set_entrez,
      array.type = "EPIC",
      sig.genes = TRUE
    )
    
    # combine with full output
    gsa_output <- gsa.region %>%
      as.data.frame() %>%
      mutate(cluster = clusters[i]) %>%
      rownames_to_column("pathway") %>%
      dplyr::arrange(FDR) %>%
      filter(FDR < 0.05)
    output_df <- rbind(output_df, gsa_output)
    
    # barplot of fgsea per cluster (plot only top 50)
    if (nrow(gsa_output) > 0) {
      gsa_output <- gsa_output %>%
        dplyr::arrange(FDR) %>%
        slice_head(n = 50)
      gsa_output$pathway <- gsub("_", " ", gsa_output$pathway)
      gsa_output$pathway <- factor(gsa_output$pathway, levels = gsa_output$pathway)
      p <- ggplot(gsa_output, aes(pathway, y = (-1) * log10(FDR), fill = FDR)) +
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
  
  # write output to tsv
  write_tsv(x = output_df, file = file.path(
    output_dir,
    paste0(prefix, "_gsaregion_output_per_cluster.tsv")
  ))
}

# run promoter analysis
run_analysis(
  methyl_m_values_full = methyl_m_values_full,
  methyl_annot_full = methyl_annot_full,
  gene_feature_filter = "promoter",
  prefix = "promoter"
)

# run gene body analysis
run_analysis(
  methyl_m_values_full = methyl_m_values_full,
  methyl_annot_full = methyl_annot_full,
  gene_feature_filter = c("exon", "intron"),
  prefix = "gene_body"
)

# run both gene body + promoter analysis
run_analysis(
  methyl_m_values_full = methyl_m_values_full,
  methyl_annot_full = methyl_annot_full,
  gene_feature_filter = c("exon", "intron", "promoter"),
  prefix = "genebody_promoter"
)
