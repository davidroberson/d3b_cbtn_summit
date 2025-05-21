#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Integrative NMF Clustering
doc: |
  Performs multi-modal clustering using the Integrative Non-negative Matrix
  Factorization (IntNMF) method on multiple data types.

baseCommand: ["Rscript", "--vanilla", "/app/cbtn_multiomics/integrative_nmf/01-multi-modal-clustering-run.R"]

requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3"
  - class: InlineJavascriptRequirement

inputs:
  samples_map:
    type: File
    inputBinding:
      position: 1
      prefix: --samples_map
  
  rna_data:
    type: File
    inputBinding:
      position: 2
      prefix: --rna_file
  
  methyl_data:
    type: File
    inputBinding:
      position: 3
      prefix: --methyl_file
  
  splice_data:
    type: File
    inputBinding:
      position: 4
      prefix: --splice_file
  
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
      prefix: --cluster_method
  
  results_dir:
    type: string
    default: "results"
    inputBinding:
      position: 7
      prefix: --results_dir
  
  plots_dir:
    type: string
    default: "plots"
    inputBinding:
      position: 8
      prefix: --plots_dir

outputs:
  cluster_assignments:
    type: File
    outputBinding:
      glob: $(inputs.results_dir)/intnmf_clusters.tsv
  
  feature_scores_rna:
    type: File
    outputBinding:
      glob: $(inputs.results_dir)/feature_scores/feature_scores_rna.tsv
  
  feature_scores_methyl:
    type: File
    outputBinding:
      glob: $(inputs.results_dir)/feature_scores/feature_scores_methyl.tsv
  
  feature_scores_splice:
    type: File
    outputBinding:
      glob: $(inputs.results_dir)/feature_scores/feature_scores_splice.tsv
  
  consensus_plot:
    type: File
    outputBinding:
      glob: $(inputs.plots_dir)/intnmf_consensus_plot.pdf
  
  silhouette_plot:
    type: File
    outputBinding:
      glob: $(inputs.plots_dir)/intnmf_silhouette_plot.pdf