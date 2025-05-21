# Function: DGE analysis by DESeq2

suppressPackageStartupMessages({
  library(tidyverse, quietly = TRUE)
  library(optparse, quietly = TRUE)
  
  # Make sure required packages are installed with clear error messages
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cran.r-project.org")
    BiocManager::install("DESeq2", update = FALSE)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cran.r-project.org")
    BiocManager::install("rtracklayer", update = FALSE)
  }
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    install.packages("msigdbr", repos="https://cran.r-project.org")
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cran.r-project.org")
    BiocManager::install("clusterProfiler", update = FALSE)
  }
  
  library(DESeq2, quietly = TRUE)
  library(rtracklayer, quietly = TRUE)
  library(msigdbr, quietly = TRUE)
  library(clusterProfiler, quietly = TRUE)
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

# Define the perform_enrichment_gsea function inline since we can't rely on working directory
perform_enrichment_gsea <- function(diffexpr_res, pathways, minGSSize, maxGSSize, prefix, plots_dir, results_dir) {
  # Create directories if they don't exist
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Check if we have any differentially expressed genes
  if (nrow(diffexpr_res) == 0) {
    cat("No differentially expressed genes found for", prefix, "\n")
    cat("No differentially expressed genes found", 
        file = file.path(results_dir, paste0(prefix, "_gsea_results.txt")))
    return(NULL)
  }
  
  # Print a sample of the gene names to help with debugging
  cat("Sample gene IDs being used:", paste(head(diffexpr_res$genes, 10), collapse=", "), "\n")
  cat("Sample pathway gene IDs:", paste(head(unique(pathways$gene), 10), collapse=", "), "\n")
  
  # Check for gene overlap between DE results and pathways
  common_genes <- intersect(diffexpr_res$genes, unique(pathways$gene))
  cat("Number of genes in common with pathway database:", length(common_genes), "\n")
  
  # If no genes overlap, create a placeholder file and exit
  if (length(common_genes) < 5) {
    cat("ERROR: Insufficient gene overlap with pathway database\n")
    cat("GSEA requires at least some overlap between gene IDs. Found only", length(common_genes), "genes in common.\n",
        "This may indicate a mismatch in gene ID formats.\n",
        "Please ensure gene IDs use the same format (e.g., gene symbols or Entrez IDs) in both datasets.\n",
        file = file.path(results_dir, paste0(prefix, "_gsea_results.txt")))
    
    # Create a placeholder output file
    placeholder_file <- file.path(results_dir, paste0(prefix, "_gsea_results.tsv"))
    cat("Pathway\tDescription\tGeneRatio\tBgRatio\tpvalue\tpadj\tqvalues\tgeneID\tCount\n", 
        file = placeholder_file)
    
    return(NULL)
  }
  
  # Prepare gene list (ranked by log2FC) - filter to only genes present in pathway database
  # to avoid "No gene can be mapped" error
  filtered_diffexpr <- diffexpr_res %>% filter(genes %in% common_genes)
  
  if (nrow(filtered_diffexpr) == 0) {
    cat("No differential genes could be mapped to pathway database\n")
    cat("No genes could be mapped to pathway database",
        file = file.path(results_dir, paste0(prefix, "_gsea_results.txt")))
    return(NULL)
  }
  
  gene_list <- filtered_diffexpr$log2FC
  names(gene_list) <- filtered_diffexpr$genes
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  # Run GSEA with error handling
  tryCatch({
    gsea_results <- clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = pathways,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = 0.05,
      eps = 0
    )
    
    # Check if results are found
    if (is.null(gsea_results) || nrow(gsea_results@result) == 0) {
      cat("No significant pathways found\n")
      # Create an empty file to indicate the analysis was run but found no results
      cat("No significant pathways found", file = file.path(results_dir, paste0(prefix, "_gsea_results.txt")))
      return(NULL)
    }
    
    # Save results
    gsea_output <- gsea_results@result
    write.table(
      gsea_output,
      file = file.path(results_dir, paste0(prefix, "_gsea_results.tsv")),
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    # Create plots (if any significant results)
    if (nrow(gsea_output) > 0) {
      # Plot top 5 upregulated and downregulated pathways
      top_pathways <- c(
        head(gsea_output[gsea_output$NES > 0, "ID"], 5),
        head(gsea_output[gsea_output$NES < 0, "ID"], 5)
      )
      
      # Create individual enrichment plots for top pathways
      for (pathway in top_pathways) {
        pdf(file.path(plots_dir, paste0(prefix, "_", pathway, "_enrichment_plot.pdf")))
        tryCatch({
          print(clusterProfiler::gseaplot2(gsea_results, geneSetID = pathway, title = pathway))
        }, error = function(e) {
          cat("Error creating plot for pathway:", pathway, "\n")
        })
        dev.off()
      }
      
      # Create summary dotplot 
      if (nrow(gsea_output) >= 3) {
        pdf(file.path(plots_dir, paste0(prefix, "_dotplot.pdf")))
        tryCatch({
          print(clusterProfiler::dotplot(gsea_results, showCategory = min(10, nrow(gsea_output))))
        }, error = function(e) {
          cat("Error creating dotplot\n")
        })
        dev.off()
      }
    }
    
    return(gsea_results)
  }, error = function(e) {
    # Handle GSEA errors gracefully
    cat("ERROR running GSEA:", conditionMessage(e), "\n")
    cat("ERROR running GSEA:", conditionMessage(e), "\n",
        file = file.path(results_dir, paste0(prefix, "_gsea_results.txt")))
    
    # Create a placeholder output file
    placeholder_file <- file.path(results_dir, paste0(prefix, "_gsea_results.tsv"))
    cat("Pathway\tDescription\tGeneRatio\tBgRatio\tpvalue\tpadj\tqvalues\tgeneID\tCount\n", 
        file = placeholder_file)
    
    return(NULL)
  })
}

# results directory
results_dir <- opt$results_dir
dir.create(results_dir, showWarnings = F, recursive = T)
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# Create subdirectories
hallmark_dir <- file.path(results_dir, "hallmark")
dir.create(hallmark_dir, showWarnings = F, recursive = T)
reactome_dir <- file.path(results_dir, "reactome")
dir.create(reactome_dir, showWarnings = F, recursive = T)

# Create placeholder files to ensure outputs exist
cat("Placeholder for pathway results", 
    file = file.path(hallmark_dir, "placeholder.tsv"))
cat("Placeholder for pathway results", 
    file = file.path(reactome_dir, "placeholder.tsv"))

# read gtf and filter to protein coding
gencode_gtf <- rtracklayer::import(con = opt$gtf_file) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  dplyr::filter(gene_type == "protein_coding") %>%
  unique()

# count data
if (grepl("\\.rds$", opt$expr_mat, ignore.case = TRUE)) {
  # For RDS files
  expr_mat <- readRDS(opt$expr_mat)
} else {
  # For TSV files
  expr_mat <- read_tsv(opt$expr_mat) %>%
    column_to_rownames(var = colnames(.)[1])
}

# read cluster information
mm_clusters <- read_tsv(file.path(opt$cluster_file))
mm_clusters <- mm_clusters %>%
  dplyr::arrange(mm_cluster) %>%
  dplyr::mutate(mm_cluster = paste0("cluster_", mm_cluster))

# get sample map information to assign sample ids
# Convert to data frame if it's a matrix (handling both cases)
if (is.matrix(expr_mat)) {
  expr_mat <- as.data.frame(expr_mat)
}

# Extract samples that match the cluster file
expr_mat <- expr_mat[, mm_clusters$Kids_First_Biospecimen_ID_RNA, drop = FALSE]
colnames(expr_mat) <- mm_clusters$sample_id

# filter expression count file to contain only protein coding gene
expr_mat <- expr_mat[rownames(expr_mat) %in% gencode_gtf$gene_name, , drop = FALSE]

# Ensure samples match annotation
expr_mat <- expr_mat[, mm_clusters$sample_id, drop = FALSE]
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
  
  # DESeq2 analysis - ensure countData is a matrix, not a data frame
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = as.matrix(round(expr_mat)),
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
  cat('Running Reactome pathway enrichment \n')
  reactome_pathways <- tryCatch({
    # Try using the new parameter names first
    msigdbr::msigdbr(collection = "C2", subcollection = "CP:REACTOME")
  }, error = function(e) {
    # Fall back to old parameter names if necessary (for older versions)
    cat("Using legacy msigdbr parameters due to:", conditionMessage(e), "\n")
    msigdbr::msigdbr(category = "C2", subcategory = "CP:REACTOME")
  })
  
  # Print diagnostic info
  cat("Number of Reactome pathways loaded:", length(unique(reactome_pathways$gs_name)), "\n")
  cat("Number of unique genes in Reactome database:", length(unique(reactome_pathways$gene_symbol)), "\n")
  
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
  cat('Running Hallmark pathway enrichment \n')
  hallmark_pathways <- tryCatch({
    # Try using the new parameter names first
    msigdbr::msigdbr(collection = "H")
  }, error = function(e) {
    # Fall back to old parameter names if necessary (for older versions)
    cat("Using legacy msigdbr parameters due to:", conditionMessage(e), "\n")
    msigdbr::msigdbr(category = "H", subcategory = NULL)
  })
  
  # Print diagnostic info
  cat("Number of Hallmark pathways loaded:", length(unique(hallmark_pathways$gs_name)), "\n")
  cat("Number of unique genes in Hallmark database:", length(unique(hallmark_pathways$gene_symbol)), "\n")
  
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
cat('Writing outputs \n')
if (nrow(output_df) > 0) {
  write_tsv(x = output_df,
            file = file.path(results_dir, "diffexpr_output_per_cluster.tsv"))
} else {
  # Create placeholder file if no DEGs found
  cat("No differentially expressed genes found\n", 
      file = file.path(results_dir, "diffexpr_output_per_cluster.tsv"))
}

# Create a plots directory to ensure it exists for CWL output
dir.create(file.path(plots_dir, "placeholder"), showWarnings = FALSE)
cat("Placeholder file", file = file.path(plots_dir, "placeholder", "placeholder.txt"))
