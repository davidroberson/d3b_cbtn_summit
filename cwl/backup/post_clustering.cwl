#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Post-Clustering Analysis
doc: |
  Performs post-clustering analysis including statistical comparisons,
  visualization, and survival analysis.

baseCommand: ["Rscript", "--vanilla"]

requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3"
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.R
        entry: |
          # Post-clustering analysis
          # Load required packages
          library(optparse)
          library(tidyverse)
          library(survival)
          library(survminer)
          library(corrplot)
          library(gridExtra)
          library(circlize)
          library(ggpubr)
          library(mclust)
          if (!require("ComplexHeatmap", quietly = TRUE))
            BiocManager::install("ComplexHeatmap")
            
          # Load libraries
          suppressPackageStartupMessages({
            library(optparse)
            library(tidyverse)
            library(survival)
            library(survminer)
            library(corrplot)
            library(ComplexHeatmap)
            library(gridExtra)
            library(circlize)
            library(ggpubr)
            library(mclust)  # For adjusted rand index
          })

          # Parse command line options
          option_list <- list(
            make_option(c("--cluster_file"), type = "character", help = "Cluster assignments (.tsv)"),
            make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
            make_option(c("--rna_data"), type = "character", help = "RNA data file (.tsv)"),
            make_option(c("--methyl_data"), type = "character", help = "Methylation data file (.tsv)"),
            make_option(c("--splice_data"), type = "character", help = "Splicing data file (.tsv)"),
            make_option(c("--feature_scores_rna"), type = "character", help = "RNA feature scores (.tsv)"),
            make_option(c("--feature_scores_methyl"), type = "character", help = "Methylation feature scores (.tsv)"),
            make_option(c("--feature_scores_splice"), type = "character", help = "Splicing feature scores (.tsv)"),
            make_option(c("--output_dir"), type = "character", help = "Path to results directory"),
            make_option(c("--plots_dir"), type = "character", help = "Path to plots directory")
          )
          opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

          # Create output directories
          output_dir <- opt$output_dir
          plots_dir <- opt$plots_dir
          dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
          dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(plots_dir, "heatmaps"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(plots_dir, "bubble_plots"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(plots_dir, "survival_plots"), showWarnings = FALSE, recursive = TRUE)
          dir.create(file.path(plots_dir, "sankey_plots"), showWarnings = FALSE, recursive = TRUE)

          # Read data files
          cat('Reading input data\n')
          clusters <- read_tsv(opt$cluster_file)
          histology <- read_tsv(opt$histology_file)
          
          # Merge histology with clusters
          sample_data <- histology %>%
            select(
              Kids_First_Biospecimen_ID, sample_id, 
              molecular_subtype, dkfz_methylation_class, sex, age_at_diagnosis,
              OS_days, OS_status, EFS_days, EFS_status, CNS_region
            ) %>%
            inner_join(clusters, by = "sample_id")
          
          # 1. Statistical associations
          cat('Performing statistical comparisons\n')
          
          # 1.1 Chi-square test of independence
          # Between multi-modal clusters and molecular subtypes
          mm_vs_subtype <- with(
            sample_data, 
            table(mm_cluster, molecular_subtype)
          )
          chisq_test <- chisq.test(mm_vs_subtype)
          
          # 1.2 Adjusted Rand Index
          ari <- adjustedRandIndex(
            sample_data$mm_cluster,
            sample_data$molecular_subtype
          )
          
          # Save statistical test results
          writeLines(
            c(
              "Chi-square test of independence between multi-modal clusters and molecular subtypes:",
              paste("X-squared =", chisq_test$statistic),
              paste("df =", chisq_test$parameter),
              paste("p-value =", chisq_test$p.value),
              "",
              "Adjusted Rand Index between multi-modal clusters and molecular subtypes:",
              paste("ARI =", ari)
            ),
            file.path(output_dir, "chisq_ari_mm_vs_subtypes.txt")
          )
          
          # 2. Create heatmaps
          cat('Creating heatmaps\n')
          
          # 2.1 Sample-level heatmaps
          # Read data matrices
          rna_data <- read_tsv(opt$rna_data) %>% column_to_rownames() %>% as.matrix()
          methyl_data <- read_tsv(opt$methyl_data) %>% column_to_rownames() %>% as.matrix()
          splice_data <- read_tsv(opt$splice_data) %>% column_to_rownames() %>% as.matrix()
          
          # Get top features per cluster from feature scores
          feature_scores_rna <- read_tsv(opt$feature_scores_rna)
          feature_scores_methyl <- read_tsv(opt$feature_scores_methyl)
          feature_scores_splice <- read_tsv(opt$feature_scores_splice)
          
          # Select top features for each data type
          top_n_features <- 50
          
          top_rna_features <- lapply(1:ncol(feature_scores_rna), function(i) {
            feature_scores_rna %>%
              select(1, !!i) %>%
              arrange(desc(!!i)) %>%
              slice_head(n = top_n_features) %>%
              pull(1)
          }) %>% unlist() %>% unique()
          
          top_methyl_features <- lapply(1:ncol(feature_scores_methyl), function(i) {
            feature_scores_methyl %>%
              select(1, !!i) %>%
              arrange(desc(!!i)) %>%
              slice_head(n = top_n_features) %>%
              pull(1)
          }) %>% unlist() %>% unique()
          
          top_splice_features <- lapply(1:ncol(feature_scores_splice), function(i) {
            feature_scores_splice %>%
              select(1, !!i) %>%
              arrange(desc(!!i)) %>%
              slice_head(n = top_n_features) %>%
              pull(1)
          }) %>% unlist() %>% unique()
          
          # Filter data matrices to top features
          rna_data_top <- rna_data[rownames(rna_data) %in% top_rna_features, ]
          methyl_data_top <- methyl_data[rownames(methyl_data) %in% top_methyl_features, ]
          splice_data_top <- splice_data[rownames(splice_data) %in% top_splice_features, ]
          
          # Create sample annotation data
          sample_annotation <- sample_data %>%
            select(sample_id, mm_cluster, molecular_subtype, dkfz_methylation_class) %>%
            column_to_rownames("sample_id")
          
          # Create sample heatmaps
          pdf(
            file = file.path(plots_dir, "heatmaps", "sample_level_heatmaps.pdf"),
            width = 12,
            height = 10
          )
          
          # Plot each data type
          # RNA heatmap
          heatmap_rna <- Heatmap(
            t(scale(t(rna_data_top))),  # Scale by feature
            name = "RNA Z-score",
            show_row_names = FALSE,
            show_column_names = FALSE,
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            column_split = sample_annotation$mm_cluster[colnames(rna_data_top)],
            column_title = "RNA Expression",
            col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
          )
          
          # Methylation heatmap
          heatmap_methyl <- Heatmap(
            t(scale(t(methyl_data_top))),  # Scale by feature
            name = "Methyl Z-score",
            show_row_names = FALSE,
            show_column_names = FALSE,
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            column_split = sample_annotation$mm_cluster[colnames(methyl_data_top)],
            column_title = "Methylation",
            col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
          )
          
          # Splicing heatmap
          heatmap_splice <- Heatmap(
            t(scale(t(splice_data_top))),  # Scale by feature
            name = "Splice Z-score",
            show_row_names = FALSE,
            show_column_names = FALSE,
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            column_split = sample_annotation$mm_cluster[colnames(splice_data_top)],
            column_title = "Splicing",
            col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
          )
          
          # Draw the heatmaps
          draw(heatmap_rna)
          draw(heatmap_methyl)
          draw(heatmap_splice)
          
          dev.off()
          
          # 2.2 Feature-level heatmaps
          # Plot feature scores
          pdf(
            file = file.path(plots_dir, "heatmaps", "feature_level_heatmaps.pdf"),
            width = 12,
            height = 10
          )
          
          # Prepare feature scores for heatmaps
          prepare_feature_scores <- function(feature_scores, top_n = 20) {
            # Get top features for each cluster
            top_features <- lapply(1:ncol(feature_scores), function(i) {
              feature_scores %>%
                select(1, !!i) %>%
                arrange(desc(!!i)) %>%
                slice_head(n = top_n) %>%
                pull(1)
            })
            
            # Combine and remove duplicates
            all_top_features <- unique(unlist(top_features))
            
            # Filter feature scores to top features
            feature_scores_mat <- as.matrix(feature_scores[
              feature_scores[[1]] %in% all_top_features, 
              -1
            ])
            rownames(feature_scores_mat) <- feature_scores[
              feature_scores[[1]] %in% all_top_features, 
              1
            ][[1]]
            
            return(feature_scores_mat)
          }
          
          # Prepare matrices
          rna_scores_mat <- prepare_feature_scores(feature_scores_rna)
          methyl_scores_mat <- prepare_feature_scores(feature_scores_methyl)
          splice_scores_mat <- prepare_feature_scores(feature_scores_splice)
          
          # Plot each feature score matrix
          heatmap_rna_scores <- Heatmap(
            rna_scores_mat,
            name = "RNA Weight",
            show_row_names = TRUE,
            show_column_names = TRUE,
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            row_names_gp = gpar(fontsize = 6),
            column_title = "RNA Feature Scores",
            col = colorRamp2(c(0, max(rna_scores_mat)/2, max(rna_scores_mat)), c("white", "pink", "red"))
          )
          
          heatmap_methyl_scores <- Heatmap(
            methyl_scores_mat,
            name = "Methyl Weight",
            show_row_names = TRUE,
            show_column_names = TRUE,
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            row_names_gp = gpar(fontsize = 6),
            column_title = "Methylation Feature Scores",
            col = colorRamp2(c(0, max(methyl_scores_mat)/2, max(methyl_scores_mat)), c("white", "pink", "red"))
          )
          
          heatmap_splice_scores <- Heatmap(
            splice_scores_mat,
            name = "Splice Weight",
            show_row_names = TRUE,
            show_column_names = TRUE,
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            row_names_gp = gpar(fontsize = 6),
            column_title = "Splicing Feature Scores",
            col = colorRamp2(c(0, max(splice_scores_mat)/2, max(splice_scores_mat)), c("white", "pink", "red"))
          )
          
          # Draw the heatmaps
          draw(heatmap_rna_scores)
          draw(heatmap_methyl_scores)
          draw(heatmap_splice_scores)
          
          dev.off()
          
          # 3. Create bubble plots
          cat('Creating bubble plots\n')
          
          # 3.1 Balloon plots
          # Multi-modal clusters vs molecular subtypes
          # Create a proper table format for data processing
          mm_vs_subtype_tab <- table(sample_data$mm_cluster, sample_data$molecular_subtype)
          mm_vs_subtype_df <- as.data.frame(mm_vs_subtype_tab)
          colnames(mm_vs_subtype_df) <- c("mm_cluster", "molecular_subtype", "Count")
          
          # Calculate percentages
          mm_vs_subtype_df <- mm_vs_subtype_df %>%
            group_by(mm_cluster) %>%
            mutate(
              Percentage = Count / sum(Count) * 100,
              Label = ifelse(Percentage > 5, sprintf("%.1f%%", Percentage), "")
            ) %>%
            ungroup()
          
          # Create bubble plot
          p_balloon <- ggplot(mm_vs_subtype_df, aes(x = mm_cluster, y = molecular_subtype)) +
            geom_point(aes(size = Percentage, color = Percentage)) +
            scale_size_continuous(range = c(1, 20)) +
            scale_color_gradient(low = "blue", high = "red") +
            geom_text(aes(label = Label), color = "white", size = 3) +
            labs(
              title = "Multi-Modal Clusters vs. Molecular Subtypes",
              x = "Multi-Modal Cluster",
              y = "Molecular Subtype",
              size = "Percentage",
              color = "Percentage"
            ) +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major = element_line(color = "gray90"),
              panel.grid.minor = element_blank()
            )
          
          # Save plot
          ggsave(
            file.path(plots_dir, "bubble_plots", "mm_clusters_vs_molsubtype_balloonplot.pdf"),
            p_balloon,
            width = 10,
            height = 8
          )
          
          # 3.2 Correlation plots
          # Multi-modal clusters vs molecular subtypes
          mm_vs_subtype_corr <- cor(
            model.matrix(~ 0 + mm_cluster, sample_data),
            model.matrix(~ 0 + molecular_subtype, sample_data)
          )
          
          pdf(
            file.path(plots_dir, "bubble_plots", "mm_clusters_vs_molsubtype_corrplot.pdf"),
            width = 10,
            height = 8
          )
          
          corrplot(
            mm_vs_subtype_corr,
            method = "circle",
            type = "full",
            tl.col = "black",
            tl.srt = 45,
            addCoef.col = "black",
            number.cex = 0.7,
            title = "Correlation between Multi-Modal Clusters and Molecular Subtypes",
            mar = c(0, 0, 2, 0)
          )
          
          dev.off()
          
          # 4. Survival analysis
          cat('Performing survival analysis\n')
          
          # 4.1 Overall survival
          # Prepare survival data
          surv_data <- sample_data %>%
            filter(!is.na(OS_days), !is.na(OS_status)) %>%
            mutate(
              OS_status_binary = ifelse(OS_status == "Dead", 1, 0),
              OS_years = OS_days / 365.25
            )
          
          # Create survival object
          surv_obj <- Surv(surv_data$OS_years, surv_data$OS_status_binary)
          
          # Fit survival curves by multi-modal cluster
          fit_mm <- survfit(surv_obj ~ mm_cluster, data = surv_data)
          
          # Perform log-rank test
          log_rank_mm <- survdiff(surv_obj ~ mm_cluster, data = surv_data)
          log_rank_p <- 1 - pchisq(log_rank_mm$chisq, df = length(unique(surv_data$mm_cluster)) - 1)
          
          # Create survival plot
          p_surv_mm <- ggsurvplot(
            fit_mm,
            data = surv_data,
            pval = TRUE,
            pval.method = TRUE,
            risk.table = TRUE,
            xlab = "Time (years)",
            legend.title = "Multi-Modal Cluster",
            title = "Overall Survival by Multi-Modal Cluster",
            ggtheme = theme_minimal()
          )
          
          # Save plot
          pdf(
            file.path(plots_dir, "survival_plots", "survival_mm_clusters.pdf"),
            width = 10,
            height = 8
          )
          print(p_surv_mm)
          dev.off()
          
          # Fit survival curves by molecular subtype
          fit_subtype <- survfit(surv_obj ~ molecular_subtype, data = surv_data)
          
          # Create survival plot
          p_surv_subtype <- ggsurvplot(
            fit_subtype,
            data = surv_data,
            pval = TRUE,
            pval.method = TRUE,
            risk.table = TRUE,
            xlab = "Time (years)",
            legend.title = "Molecular Subtype",
            title = "Overall Survival by Molecular Subtype",
            ggtheme = theme_minimal()
          )
          
          # Save plot
          pdf(
            file.path(plots_dir, "survival_plots", "survival_molsubtype.pdf"),
            width = 10,
            height = 8
          )
          print(p_surv_subtype)
          dev.off()
          
          # 4.2 Cox proportional hazards model
          # Fit Cox model
          cox_mm <- coxph(surv_obj ~ mm_cluster, data = surv_data)
          cox_summary <- summary(cox_mm)
          
          # Save model summary
          capture.output(
            cox_summary,
            file = file.path(output_dir, "coxph_summary_OS.txt")
          )
          
          # Calculate risk scores
          risk_scores <- data.frame(
            sample_id = surv_data$sample_id,
            mm_cluster = surv_data$mm_cluster,
            risk_score = predict(cox_mm)
          )
          
          # Save risk scores
          write_tsv(
            risk_scores,
            file = file.path(output_dir, "coxph_risk_score_OS.txt")
          )
          
          # Plot risk scores
          p_risk <- ggplot(risk_scores, aes(x = mm_cluster, y = risk_score, fill = mm_cluster)) +
            geom_boxplot() +
            labs(
              title = "Risk Scores by Multi-Modal Cluster",
              x = "Multi-Modal Cluster",
              y = "Risk Score (log hazard ratio)"
            ) +
            theme_minimal()
          
          # Save plot
          pdf(
            file.path(plots_dir, "survival_plots", "coxph_summary_OS.pdf"),
            width = 10,
            height = 8
          )
          print(p_risk)
          dev.off()
          
          cat('Post-clustering analysis complete\n')

inputs:
  cluster_file:
    type: File
    inputBinding:
      position: 1
      prefix: --cluster_file
  
  histology_file:
    type: File
    inputBinding:
      position: 2
      prefix: --histology_file
  
  rna_data:
    type: File
    inputBinding:
      position: 3
      prefix: --rna_data
  
  methyl_data:
    type: File
    inputBinding:
      position: 4
      prefix: --methyl_data
  
  splice_data:
    type: File
    inputBinding:
      position: 5
      prefix: --splice_data
  
  feature_scores_rna:
    type: File
    inputBinding:
      position: 6
      prefix: --feature_scores_rna
  
  feature_scores_methyl:
    type: File
    inputBinding:
      position: 7
      prefix: --feature_scores_methyl
  
  feature_scores_splice:
    type: File
    inputBinding:
      position: 8
      prefix: --feature_scores_splice
  
  output_dir:
    type: string
    default: "results"
    inputBinding:
      position: 9
      prefix: --output_dir
  
  plots_dir:
    type: string
    default: "plots"
    inputBinding:
      position: 10
      prefix: --plots_dir

arguments:
  - position: 0
    valueFrom: "script.R"

outputs:
  comparison_results:
    type: File
    outputBinding:
      glob: "$(inputs.output_dir)/chisq_ari_mm_vs_subtypes.txt"
  
  heatmaps:
    type: Directory
    outputBinding:
      glob: "$(inputs.plots_dir)/heatmaps"
  
  bubble_plots:
    type: Directory
    outputBinding:
      glob: "$(inputs.plots_dir)/bubble_plots"
  
  survival_plots:
    type: Directory
    outputBinding:
      glob: "$(inputs.plots_dir)/survival_plots"
  
  sankey_plots:
    type: Directory
    outputBinding:
      glob: "$(inputs.plots_dir)/sankey_plots"