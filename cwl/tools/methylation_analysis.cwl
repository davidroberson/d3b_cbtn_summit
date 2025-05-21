#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Methylation Analysis
doc: |
  Performs differential methylation analysis between clusters
  and conducts pathway enrichment analysis on the differentially
  methylated sites.

baseCommand: ["Rscript", "--vanilla", "/app/cbtn_multiomics/methylation/01-limma_analysis.R"]

requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3"
  - class: InlineJavascriptRequirement

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
  
  results_dir:
    type: string
    default: "results"
    inputBinding:
      position: 3
      prefix: --results_dir
  
  plots_dir:
    type: string
    default: "plots"
    inputBinding:
      position: 4
      prefix: --plots_dir

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