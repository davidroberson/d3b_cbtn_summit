#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Differential Gene Expression Analysis
doc: |
  Performs differential gene expression analysis and pathway enrichment
  for each cluster versus the rest.

baseCommand: ["Rscript", "--vanilla"]

requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.1"
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.R
        entry: |
          # Differential gene expression analysis
          # Load required packages
          library(optparse)
          library(tidyverse)
          library(ggplot2)
          library(ggpubr)
          library(msigdbr)
          library(DESeq2)
          library(rtracklayer)
          library(clusterProfiler)
          
          # Suppress package startup messages
          suppressPackageStartupMessages({
            library(optparse)
            library(tidyverse)
            library(DESeq2)
            library(rtracklayer)
            library(clusterProfiler)
            library(msigdbr)
            library(ggplot2)
            library(ggpubr)
          })

          # Parse command line options
          option_list <- list(
            make_option(c("--expr_mat"), type = "character", help = "Expression matrix, counts (.rds)"),
            make_option(c("--gtf_file"), type = "character", help = "Gencode gtf file"),
            make_option(c("--cluster_file"), type = "character", help = "Cluster assignments (.tsv)"),
            make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
            make_option(c("--results_dir"), type = "character", help = "Path to results directory"),
            make_option(c("--plots_dir"), type = "character", help = "Path to plots directory")
          )
          opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

          # Create output directories
          results_dir <- opt$results_dir
          plots_dir <- opt$plots_dir
          dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
          dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(results_dir, "hallmark"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(results_dir, "reactome"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(plots_dir, "hallmark"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(plots_dir, "reactome"), showWarnings = FALSE, recursive = TRUE)

          # Read gene annotation
          cat('Reading gene annotation\n')
          gtf_file <- opt$gtf_file
          gencode_gtf <- rtracklayer::import(con = gtf_file) %>%
            as.data.frame() %>%
            dplyr::select(gene_id, gene_name, gene_type) %>%
            dplyr::filter(!grepl("ENSG", gene_name), gene_type == "protein_coding") %>%
            unique()

          # Read histology file
          cat('Reading histology file\n')
          histology_file <- read_tsv(opt$histology_file)

          # Read cluster assignments
          cat('Reading cluster assignments\n')
          cluster_data <- read_tsv(opt$cluster_file)
          
          # Get number of clusters
          clusters <- sort(unique(cluster_data$mm_cluster))
          
          # Read expression data
          cat('Reading expression data\n')
          expr_mat <- readRDS(opt$expr_mat)
          
          # Filter to match samples in cluster data
          expr_mat <- expr_mat[, colnames(expr_mat) %in% histology_file$Kids_First_Biospecimen_ID]
          
          # Create a sample mapping
          sample_map <- histology_file %>%
            dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
            dplyr::filter(sample_id %in% cluster_data$sample_id)
          
          # Create sample metadata for DESeq2
          sample_metadata <- sample_map %>%
            inner_join(cluster_data, by = "sample_id")
          
          # Initialize results data frame
          all_results <- data.frame()
          
          # For each cluster, perform differential expression analysis
          for (cluster in clusters) {
            cat(sprintf('Performing differential expression for cluster %s\n', cluster))
            
            # Create comparison: this cluster vs all others
            sample_metadata$is_cluster <- ifelse(sample_metadata$mm_cluster == cluster, 1, 0)
            sample_metadata$is_cluster <- factor(sample_metadata$is_cluster)
            
            # Create DESeq2 object
            dds <- DESeqDataSetFromMatrix(
              countData = expr_mat[, sample_metadata$Kids_First_Biospecimen_ID],
              colData = sample_metadata,
              design = ~ is_cluster
            )
            
            # Run DESeq2
            dds <- DESeq(dds)
            
            # Get results
            res <- results(dds, contrast = c("is_cluster", "1", "0"))
            res_df <- as.data.frame(res) %>%
              rownames_to_column("gene_name") %>%
              dplyr::filter(!is.na(padj))
            
            # Add cluster information
            res_df$cluster <- cluster
            
            # Append to all results
            all_results <- bind_rows(all_results, res_df)
            
            # Filter for significant genes
            sig_genes <- res_df %>%
              dplyr::filter(padj < 0.05) %>%
              arrange(padj)
            
            # Perform GSEA on all genes, ranked by log2FoldChange
            ranked_genes <- res_df$log2FoldChange
            names(ranked_genes) <- res_df$gene_name
            ranked_genes <- sort(ranked_genes, decreasing = TRUE)
            
            # Get pathway gene sets
            hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")
            reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
            
            # Run GSEA for hallmark
            gsea_hallmark <- GSEA(
              geneList = ranked_genes,
              TERM2GENE = hallmark_gene_sets %>% dplyr::select(gs_name, gene_symbol),
              pvalueCutoff = 0.1,
              pAdjustMethod = "BH",
              verbose = FALSE
            )
            
            # Run GSEA for reactome
            gsea_reactome <- GSEA(
              geneList = ranked_genes,
              TERM2GENE = reactome_gene_sets %>% dplyr::select(gs_name, gene_symbol),
              pvalueCutoff = 0.1,
              pAdjustMethod = "BH",
              verbose = FALSE
            )
            
            # Save results
            if (!is.null(gsea_hallmark) && nrow(gsea_hallmark@result) > 0) {
              write_tsv(
                gsea_hallmark@result,
                file = file.path(results_dir, "hallmark", sprintf("cluster_%s_vs_rest_gsea.tsv", cluster))
              )
              
              # Create plots
              # Barplot
              pdf(
                file = file.path(plots_dir, "hallmark", sprintf("cluster_%s_vs_rest_gsea_barplot.pdf", cluster)),
                width = 10,
                height = 8
              )
              print(barplot(gsea_hallmark, showCategory = 20))
              dev.off()
              
              # Dotplot
              pdf(
                file = file.path(plots_dir, "hallmark", sprintf("cluster_%s_vs_rest_gsea_dotplot.pdf", cluster)),
                width = 10,
                height = 8
              )
              print(dotplot(gsea_hallmark, showCategory = 20))
              dev.off()
            }
            
            if (!is.null(gsea_reactome) && nrow(gsea_reactome@result) > 0) {
              write_tsv(
                gsea_reactome@result,
                file = file.path(results_dir, "reactome", sprintf("cluster_%s_vs_rest_gsea.tsv", cluster))
              )
              
              # Create plots
              # Barplot
              pdf(
                file = file.path(plots_dir, "reactome", sprintf("cluster_%s_vs_rest_gsea_barplot.pdf", cluster)),
                width = 10,
                height = 8
              )
              print(barplot(gsea_reactome, showCategory = 20))
              dev.off()
              
              # Dotplot
              pdf(
                file = file.path(plots_dir, "reactome", sprintf("cluster_%s_vs_rest_gsea_dotplot.pdf", cluster)),
                width = 10,
                height = 8
              )
              print(dotplot(gsea_reactome, showCategory = 20))
              dev.off()
            }
          }

          # Save all differential expression results
          write_tsv(
            all_results,
            file = file.path(results_dir, "diffexpr_output_per_cluster.tsv")
          )

          cat('Differential gene expression analysis complete\n')

inputs:
  expr_mat:
    type: File
    inputBinding:
      position: 1
      prefix: --expr_mat
  
  gtf_file:
    type: File
    inputBinding:
      position: 2
      prefix: --gtf_file
  
  cluster_file:
    type: File
    inputBinding:
      position: 3
      prefix: --cluster_file
  
  histology_file:
    type: File
    inputBinding:
      position: 4
      prefix: --histology_file
  
  results_dir:
    type: string
    default: "results"
    inputBinding:
      position: 5
      prefix: --results_dir
  
  plots_dir:
    type: string
    default: "plots"
    inputBinding:
      position: 6
      prefix: --plots_dir

arguments:
  - position: 0
    valueFrom: "script.R"

outputs:
  deseq_results:
    type: File
    outputBinding:
      glob: $(inputs.results_dir)/diffexpr_output_per_cluster.tsv
  
  pathway_results:
    type: File[]
    outputBinding:
      glob: [
        "$(inputs.results_dir)/hallmark/*.tsv",
        "$(inputs.results_dir)/reactome/*.tsv"
      ]
  
  pathway_plots:
    type: Directory
    outputBinding:
      glob: $(inputs.plots_dir)