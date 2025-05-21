#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Differential Gene Expression Analysis
doc: |
  Performs differential gene expression analysis between clusters
  and the rest of the samples using DESeq2, and conducts pathway
  enrichment analysis on the results.

baseCommand: ["Rscript", "--vanilla", "/app/cbtn_multiomics/dge_analysis/01-dge_analysis_deseq.R"]

requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3"
  - class: InlineJavascriptRequirement

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
        "$(inputs.results_dir)/reactome/*.tsv",
        "*.tsv_*"
      ]
  
  pathway_plots:
    type: Directory
    outputBinding:
      glob: $(inputs.plots_dir)