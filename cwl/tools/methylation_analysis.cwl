#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Methylation Analysis
doc: |
  Performs differential methylation analysis and pathway enrichment
  for each cluster versus the rest.

baseCommand: ["Rscript", "--vanilla"]

requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0"
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.R
        entry: |
          # Methylation analysis
          # Load required packages
          suppressPackageStartupMessages({
            library(optparse)
            library(tidyverse)
            library(limma)
            library(missMethyl)
            library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
            library(ggplot2)
          })

          # Parse command line options
          option_list <- list(
            make_option(c("--methyl_file"), type = "character", help = "Methylation m-values (.rds)"),
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
          dir.create(file.path(results_dir, "limma_output"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(results_dir, "dms_gsameth_output"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(results_dir, "dms_gsameth_output", "hallmark"), showWarnings = FALSE, recursive = TRUE)
          dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(plots_dir, "dms_gsameth_output"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(plots_dir, "dms_gsameth_output", "hallmark"), showWarnings = FALSE, recursive = TRUE)

          # Read histology file
          cat('Reading histology file\n')
          histology_file <- read_tsv(opt$histology_file)

          # Read cluster assignments
          cat('Reading cluster assignments\n')
          cluster_data <- read_tsv(opt$cluster_file)
          
          # Get number of clusters
          clusters <- sort(unique(cluster_data$mm_cluster))
          
          # Read methylation data
          cat('Reading methylation data\n')
          methyl_data <- readRDS(opt$methyl_file)
          
          # Get EPIC annotation
          anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
          
          # Filter for probes in promoter regions
          promoter_probes <- anno$Name[anno$UCSC_RefGene_Group %in% c("TSS1500", "TSS200", "1stExon")]
          genebody_probes <- anno$Name[anno$UCSC_RefGene_Group %in% c("Body", "3'UTR", "5'UTR")]
          combined_probes <- anno$Name[anno$UCSC_RefGene_Group %in% c("TSS1500", "TSS200", "1stExon", "Body", "3'UTR", "5'UTR")]
          
          # Filter methylation data
          methyl_promoter <- methyl_data[rownames(methyl_data) %in% promoter_probes, ]
          methyl_genebody <- methyl_data[rownames(methyl_data) %in% genebody_probes, ]
          methyl_combined <- methyl_data[rownames(methyl_data) %in% combined_probes, ]
          
          # Create a mapping between samples and clusters
          sample_map <- histology_file %>%
            dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
            dplyr::filter(sample_id %in% cluster_data$sample_id) %>%
            inner_join(cluster_data, by = "sample_id")
          
          # Filter methylation data to match samples
          methyl_promoter <- methyl_promoter[, colnames(methyl_promoter) %in% sample_map$Kids_First_Biospecimen_ID]
          methyl_genebody <- methyl_genebody[, colnames(methyl_genebody) %in% sample_map$Kids_First_Biospecimen_ID]
          methyl_combined <- methyl_combined[, colnames(methyl_combined) %in% sample_map$Kids_First_Biospecimen_ID]
          
          # Initialize results data frames
          promoter_results <- data.frame()
          genebody_results <- data.frame()
          combined_results <- data.frame()
          
          # For each cluster, perform differential methylation analysis
          for (cluster in clusters) {
            cat(sprintf('Performing differential methylation for cluster %s\n', cluster))
            
            # Create design matrix
            is_cluster <- ifelse(
              sample_map$mm_cluster[match(colnames(methyl_promoter), sample_map$Kids_First_Biospecimen_ID)] == cluster,
              1, 0
            )
            design <- model.matrix(~ is_cluster)
            
            # Perform limma analysis for promoter probes
            fit_promoter <- lmFit(methyl_promoter, design)
            fit_promoter <- eBayes(fit_promoter)
            tt_promoter <- topTable(fit_promoter, coef = 2, number = Inf)
            tt_promoter$cluster <- cluster
            promoter_results <- bind_rows(promoter_results, tt_promoter)
            
            # Perform limma analysis for gene body probes
            fit_genebody <- lmFit(methyl_genebody, design)
            fit_genebody <- eBayes(fit_genebody)
            tt_genebody <- topTable(fit_genebody, coef = 2, number = Inf)
            tt_genebody$cluster <- cluster
            genebody_results <- bind_rows(genebody_results, tt_genebody)
            
            # Perform limma analysis for combined probes
            fit_combined <- lmFit(methyl_combined, design)
            fit_combined <- eBayes(fit_combined)
            tt_combined <- topTable(fit_combined, coef = 2, number = Inf)
            tt_combined$cluster <- cluster
            combined_results <- bind_rows(combined_results, tt_combined)
            
            # For combined analysis, perform gsameth for pathway analysis
            sig_dmps <- rownames(tt_combined)[tt_combined$adj.P.Val < 0.05]
            if (length(sig_dmps) > 10) {
              # If too many DMPs, take top 10000
              if (length(sig_dmps) > 10000) {
                sig_dmps <- rownames(tt_combined)[order(tt_combined$adj.P.Val)][1:10000]
              }
              
              # Get universe of probes
              all_probes <- rownames(methyl_combined)
              
              # Perform gsameth for hallmark gene sets
              gsameth_hall <- gsameth(
                sig.cpg = sig_dmps,
                all.cpg = all_probes,
                collection = "h.all.v2023.2",
                array.type = "EPIC"
              )
              
              # Save results
              if (!is.null(gsameth_hall) && nrow(gsameth_hall) > 0) {
                write_tsv(
                  gsameth_hall,
                  file = file.path(
                    results_dir, "dms_gsameth_output", "hallmark",
                    sprintf("cluster_%s_vs_rest_gsameth.tsv", cluster)
                  )
                )
                
                # Create bar plot for top pathways
                if (nrow(gsameth_hall) > 0) {
                  pdf(
                    file = file.path(
                      plots_dir, "dms_gsameth_output", "hallmark",
                      sprintf("cluster_%s_vs_rest_gsameth_barplot.pdf", cluster)
                    ),
                    width = 10,
                    height = 8
                  )
                  top_pathways <- gsameth_hall %>%
                    filter(FDR < 0.1) %>%
                    arrange(FDR) %>%
                    head(50)
                  if (nrow(top_pathways) > 0) {
                    p <- ggplot(top_pathways, aes(x = reorder(ID, -log10(FDR)), y = -log10(FDR))) +
                      geom_bar(stat = "identity", fill = "steelblue") +
                      coord_flip() +
                      labs(
                        title = sprintf("Cluster %s vs Rest: GSA Methylation (Hallmark)", cluster),
                        x = "Pathway",
                        y = "-log10(FDR)"
                      ) +
                      theme_minimal() +
                      theme(axis.text.y = element_text(size = 8))
                    print(p)
                  }
                  dev.off()
                }
              }
            }
          }

          # Save all differential methylation results
          write_tsv(
            promoter_results,
            file = file.path(results_dir, "limma_output", "promoter_diffexpr_probes_per_cluster.tsv")
          )
          
          write_tsv(
            genebody_results,
            file = file.path(results_dir, "limma_output", "gene_body_diffexpr_probes_per_cluster.tsv")
          )
          
          write_tsv(
            combined_results,
            file = file.path(results_dir, "limma_output", "genebody_promoter_diffexpr_probes_per_cluster.tsv")
          )
          
          # Combine all gsameth results
          gsameth_files <- list.files(
            path = file.path(results_dir, "dms_gsameth_output", "hallmark"),
            pattern = ".*_gsameth.tsv$",
            full.names = TRUE
          )
          
          if (length(gsameth_files) > 0) {
            all_gsameth <- lapply(gsameth_files, function(file) {
              cluster <- gsub(".*cluster_(.*)_vs_rest_gsameth.tsv", "\\1", basename(file))
              result <- read_tsv(file)
              result$cluster <- cluster
              return(result)
            }) %>% bind_rows()
            
            write_tsv(
              all_gsameth,
              file = file.path(
                results_dir, "dms_gsameth_output", "hallmark",
                "genebody_promoter_gsameth_output_per_cluster.tsv"
              )
            )
            
            # Create a summary plot
            pdf(
              file = file.path(
                plots_dir, "dms_gsameth_output", "hallmark",
                "genebody_promoter_gsameth_pathways.pdf"
              ),
              width = 12,
              height = 10
            )
            
            top_pathways_per_cluster <- all_gsameth %>%
              filter(FDR < 0.1) %>%
              group_by(cluster) %>%
              arrange(FDR) %>%
              slice_head(n = 10) %>%
              ungroup()
            
            if (nrow(top_pathways_per_cluster) > 0) {
              p <- ggplot(top_pathways_per_cluster, aes(x = reorder(ID, -log10(FDR)), y = -log10(FDR), fill = cluster)) +
                geom_bar(stat = "identity", position = "dodge") +
                coord_flip() +
                labs(
                  title = "Top Methylation Pathways by Cluster",
                  x = "Pathway",
                  y = "-log10(FDR)",
                  fill = "Cluster"
                ) +
                theme_minimal() +
                theme(axis.text.y = element_text(size = 8))
              print(p)
            }
            
            dev.off()
          }

          cat('Methylation analysis complete\n')

inputs:
  methyl_file:
    type: File
    inputBinding:
      position: 1
      prefix: --methyl_file
  
  cluster_file:
    type: File
    inputBinding:
      position: 2
      prefix: --cluster_file
  
  histology_file:
    type: File
    inputBinding:
      position: 3
      prefix: --histology_file
  
  results_dir:
    type: string
    default: "results"
    inputBinding:
      position: 4
      prefix: --results_dir
  
  plots_dir:
    type: string
    default: "plots"
    inputBinding:
      position: 5
      prefix: --plots_dir

arguments:
  - position: 0
    valueFrom: "script.R"

outputs:
  methylation_results:
    type: File[]
    outputBinding:
      glob: "$(inputs.results_dir)/limma_output/*.tsv"
  
  pathway_results:
    type: File
    outputBinding:
      glob: "$(inputs.results_dir)/dms_gsameth_output/hallmark/genebody_promoter_gsameth_output_per_cluster.tsv"
  
  pathway_plots:
    type: Directory
    outputBinding:
      glob: $(inputs.plots_dir)