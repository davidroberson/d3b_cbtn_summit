#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Integrative NMF Clustering
doc: |
  Performs integrative non-negative matrix factorization (IntNMF)
  on multiple data modalities to discover clusters.

baseCommand: ["Rscript", "--vanilla"]

requirements:
  - class: DockerRequirement
    # Use our custom image with all required packages pre-installed
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0"
    # When testing with sbpack, use your own repository
    # dockerPull: "<your-dockerhub-username>/cbtn-summit-workflow:1.0.0"
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.R
        entry: |
          #!/usr/bin/env Rscript
          # run multi-modal clustering
          
          # Load required packages
          library(optparse)
          library(tidyverse)
          library(IntNMF)
          library(matrixStats)
          library(cluster)
          
          suppressPackageStartupMessages({
            library(optparse)
            library(IntNMF)
            library(tidyverse)
            library(matrixStats)
            library(cluster)
          })
          
          # parse command line options
          option_list <- list(
            make_option(c("--rna_data"), type = "character", help = "RNA data file (.tsv)"),
            make_option(c("--methyl_data"), type = "character", help = "Methylation data file (.tsv)"),
            make_option(c("--splice_data"), type = "character", help = "Splicing data file (.tsv)"),
            make_option(c("--samples_map"), type = "character", help = "Sample mapping file (.tsv)"),
            make_option(c("--output_dir"), type = "character", help = "Path to results directory"),
            make_option(c("--plots_dir"), type = "character", help = "Path to plots directory"),
            make_option(c("--max_k"), type = "integer", default = 10, help = "Maximum number of clusters to try"),
            make_option(c("--method"), type = "character", default = "intNMF", help = "Clustering method to use")
          )
          opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
          max_k <- opt$max_k
          method <- opt$method
          
          # Log which method is being used
          cat('Using clustering method:', method, '\n')
          
          # output directory
          output_dir <- opt$output_dir
          dir.create(output_dir, showWarnings = F, recursive = T)
          
          # plots directory
          plots_dir <- opt$plots_dir
          dir.create(plots_dir, showWarnings = F, recursive = T)
          
          # source function
          source("utils/run_clusterstats.R")
          
          # read data matrices and add to a list
          dat <- list()
          
          count_data <- read_tsv(opt$rna_data) %>%
            column_to_rownames() %>%
            as.matrix()
          dat[["rna_data"]] <- count_data
          methyl_data <- read_tsv(opt$methyl_data) %>%
            column_to_rownames() %>%
            as.matrix()
          dat[["methyl_data"]] <- methyl_data
          splice_data <- read_tsv(opt$splice_data) %>%
            column_to_rownames() %>%
            as.matrix()
          dat[["splice_data"]] <- splice_data
          
          # add weights to each data modality
          wt = if (is.list(dat))
            rep(1, length(dat)) else 1
          
          # get the nmf output corresponding to the most optimal k
          cat('Running multi-modal clustering \n')
          nmf_output <- run_clusterstats(
            dat = dat,
            output_dir = output_dir,
            k_value = max_k
          )
          
          # save feature scores per cluster for downstream processing
          subdir <- file.path(output_dir, "feature_scores")
          dir.create(subdir, showWarnings = F, recursive = T)
          write_tsv(as.data.frame(nmf_output$H$rna_data),
                    file = file.path(subdir, "feature_scores_rna.tsv"))
          write_tsv(
            as.data.frame(nmf_output$H$methyl_data),
            file = file.path(subdir, "feature_scores_methyl.tsv")
          )
          write_tsv(
            as.data.frame(nmf_output$H$splice_data),
            file = file.path(subdir, "feature_scores_splice.tsv")
          )
          
          # 1) ConsensusMatPlot
          # Given the integrative NMF fit object, the function creates image plot of the consensus matrix ordered
          # according to clusters groups. Cleaner block structure indicates stronger clusters
          cat('Generating cluster quality plots \n')
          pdf(
            file = file.path(plots_dir, "intnmf_consensus_plot.pdf"),
            width = 10,
            height = 10,
            onefile = F
          )
          ConsensusMatPlot(fit = nmf_output,
                         rowLab = TRUE,
                         colLab = TRUE)
          dev.off()
          
          # 2) SilhouettePlot
          # Silhouette width plot is returned together with mean silhouette width for each group, overall silhouette width and summary statistics.
          pdf(
            file = file.path(plots_dir, "intnmf_silhouette_plot.pdf"),
            width = 10,
            height = 10,
            onefile = F
          )
          SilhouettePlot(fit = nmf_output, cluster.col = NULL)
          dev.off()
          
          # 3) output clusters with per sample
          df <- data.frame(sample_id = names(nmf_output$clusters),
                         mm_cluster = nmf_output$clusters)
          
          # combine with sample mapping file that has biospecimens for each data modality
          samples_map <- read_tsv(opt$samples_map)
          df <- samples_map %>%
            inner_join(df)
          
          # write output
          cat('Writing clustering output \n')
          write_tsv(df, file = file.path(output_dir, "intnmf_clusters.tsv"))
      - entryname: utils/run_clusterstats.R
        entry: |
          # This script is to calculate cluster stats for each k
          suppressPackageStartupMessages({
            library(IntNMF)
            library(tidyverse)
          })
          
          run_clusterstats <- function(dat, output_dir, k_value) {
            # if weight is not assigned, use default
            if (is.null(wt)) {
              wt = if (is.list(dat))
                rep(1, length(dat))
              else
                1
            }
            
            # do this for each cluster
            # run Nonnegative Matrix Factorization of Multiple data using Nonnegative Alternating Least Square
            nmf_fname <- file.path(output_dir, "intnmf_fit_all.rds")
            
            
            # Normalize by each omics type's frobenius norm
            count_data_norm <- dat[["rna_data"]] / norm(as.matrix(dat[["rna_data"]]), type="F")
            methyl_data_norm <- dat[["methyl_data"]] / norm(as.matrix(dat[["methyl_data"]]), type="F")
            splice_data_norm <- dat[["splice_data"]] / norm(as.matrix(dat[["splice_data"]]), type="F")
            
            nmf_output <- nmf.mnnals(dat = list(count_data_norm, methyl_data_norm, splice_data_norm), 
                                   k = k_value, 
                                   maxiter = 200, 
                                   st.count = 20, 
                                   n.ini = 30, 
                                   ini.nndsvd = TRUE, 
                                   seed = TRUE,
                                   wt=if(is.list(dat)) rep(1,length(dat)) else 1)
            
            # add modality names instead of H1, H2... for better clarity in downstream plotting
            names(nmf_output$H) <- names(dat)
            
            # save best fit to file
            saveRDS(object = nmf_output,
                    file = file.path(output_dir, "intnmf_best_fit.rds"))
            
            # return the nmf output corresponding to the most optimal k for downstream analyses
            return(nmf_output)
          }

inputs:
  rna_data:
    type: File
    inputBinding:
      position: 1
      prefix: --rna_data
  
  methyl_data:
    type: File
    inputBinding:
      position: 2
      prefix: --methyl_data
  
  splice_data:
    type: File
    inputBinding:
      position: 3
      prefix: --splice_data
  
  samples_map:
    type: File
    inputBinding:
      position: 4
      prefix: --samples_map
  
  max_k:
    type: int
    default: 10
    inputBinding:
      position: 5
      prefix: --max_k
      
  method:
    type: string
    default: "intNMF"
    inputBinding:
      position: 6
      prefix: --method
  
  output_dir:
    type: string
    default: "results"
    inputBinding:
      position: 7
      prefix: --output_dir
  
  plots_dir:
    type: string
    default: "plots"
    inputBinding:
      position: 8
      prefix: --plots_dir

arguments:
  - position: 0
    valueFrom: "script.R"

outputs:
  cluster_assignments:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/intnmf_clusters.tsv
  
  feature_scores_rna:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/feature_scores/feature_scores_rna.tsv
  
  feature_scores_methyl:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/feature_scores/feature_scores_methyl.tsv
  
  feature_scores_splice:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/feature_scores/feature_scores_splice.tsv
  
  consensus_plot:
    type: File
    outputBinding:
      glob: $(inputs.plots_dir)/intnmf_consensus_plot.pdf
  
  silhouette_plot:
    type: File
    outputBinding:
      glob: $(inputs.plots_dir)/intnmf_silhouette_plot.pdf