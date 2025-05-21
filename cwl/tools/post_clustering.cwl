#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Post-Clustering Analysis
doc: |
  Performs post-clustering analysis including statistical comparisons,
  visualization, and survival analysis.

baseCommand: ["Rscript", "--vanilla", "/app/cbtn_multiomics/post_clustering/01-compare-classes.R"]

requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3"
  - class: InlineJavascriptRequirement

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
      prefix: --rna_file
  
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