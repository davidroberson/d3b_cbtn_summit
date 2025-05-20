#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: Workflow

label: CBTN Multi-Omic Clustering Workflow
doc: |
  A workflow for multi-modal clustering analysis of cancer genomics data.
  Takes RNA-seq expression, methylation, and alternative splicing data
  as input to generate integrated clusters and perform downstream analyses.

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement

inputs:
  histology_file:
    type: File
    doc: "Clinical and histology data for samples (.tsv)"
  
  short_histology:
    type: string
    doc: "Short histology code to filter samples (e.g., 'HGAT')"
  
  count_file:
    type: File
    doc: "RNA-seq expression matrix (.rds)"
  
  methyl_file:
    type: File
    doc: "Methylation beta values matrix (.rds)"
  
  splice_file:
    type: File
    doc: "Alternative splice junction PSI values (.rds)"
  
  gtf_file:
    type: File
    doc: "Gene annotation GTF file (.gtf.gz)"
  
  # Optional parameters with defaults
  cluster_method:
    type: string
    default: "intNMF"
    doc: "Clustering method to use (currently only intNMF supported)"
  
  num_features:
    type: int
    default: 1000
    doc: "Number of top variable features to select for each data modality"
  
  max_k:
    type: int
    default: 10
    doc: "Maximum number of clusters to try"

steps:
  data_preparation:
    run: tools/data_preparation.cwl
    in:
      histology_file: histology_file
      short_histology: short_histology
      count_file: count_file
      methyl_file: methyl_file
      splice_file: splice_file
      gtf_file: gtf_file
      num_features: num_features
    out:
      - rna_data
      - methyl_data
      - splice_data
      - samples_map

  cluster_analysis:
    run: tools/integrative_nmf.cwl
    in:
      samples_map: data_preparation/samples_map
      rna_data: data_preparation/rna_data
      methyl_data: data_preparation/methyl_data
      splice_data: data_preparation/splice_data
      max_k: max_k
      method: cluster_method
    out:
      - cluster_assignments
      - feature_scores_rna
      - feature_scores_methyl
      - feature_scores_splice
      - consensus_plot
      - silhouette_plot

  differential_expression:
    run: tools/dge_analysis.cwl
    in:
      expr_mat: count_file
      gtf_file: gtf_file
      cluster_file: cluster_analysis/cluster_assignments
      histology_file: histology_file
    out:
      - deseq_results
      - pathway_results
      - pathway_plots

  methylation_analysis:
    run: tools/methylation_analysis.cwl
    in:
      methyl_file: methyl_file
      cluster_file: cluster_analysis/cluster_assignments
      histology_file: histology_file
    out:
      - methylation_results
      - pathway_results
      - pathway_plots

  post_clustering:
    run: tools/post_clustering.cwl
    in:
      cluster_file: cluster_analysis/cluster_assignments
      histology_file: histology_file
      rna_data: data_preparation/rna_data
      methyl_data: data_preparation/methyl_data
      splice_data: data_preparation/splice_data
      feature_scores_rna: cluster_analysis/feature_scores_rna
      feature_scores_methyl: cluster_analysis/feature_scores_methyl
      feature_scores_splice: cluster_analysis/feature_scores_splice
    out:
      - comparison_results
      - heatmaps
      - bubble_plots
      - survival_plots
      - sankey_plots

outputs:
  # Organized outputs with logical grouping
  
  # 1. Processed Data Files
  processed_rna_data:
    type: File
    doc: "Processed RNA-seq data matrix"
    outputSource: data_preparation/rna_data
    
  processed_methyl_data:
    type: File
    doc: "Processed methylation data matrix"
    outputSource: data_preparation/methyl_data
    
  processed_splice_data:
    type: File
    doc: "Processed splicing data matrix"
    outputSource: data_preparation/splice_data
    
  sample_mapping:
    type: File
    doc: "Sample ID mapping file"
    outputSource: data_preparation/samples_map
  
  # 2. Clustering Results
  cluster_assignments:
    type: File
    doc: "Cluster assignments for each sample"
    outputSource: cluster_analysis/cluster_assignments
  
  feature_scores_rna:
    type: File
    doc: "RNA feature importance scores"
    outputSource: cluster_analysis/feature_scores_rna
    
  feature_scores_methyl:
    type: File
    doc: "Methylation feature importance scores"
    outputSource: cluster_analysis/feature_scores_methyl
    
  feature_scores_splice:
    type: File
    doc: "Splicing feature importance scores"
    outputSource: cluster_analysis/feature_scores_splice
  
  consensus_plot:
    type: File
    doc: "Consensus matrix plot for clustering"
    outputSource: cluster_analysis/consensus_plot
    
  silhouette_plot:
    type: File
    doc: "Silhouette plot for clustering"
    outputSource: cluster_analysis/silhouette_plot
  
  # 3. Differential Expression Results
  deseq_results:
    type: File
    doc: "Differential expression analysis results"
    outputSource: differential_expression/deseq_results
  
  expression_pathway_results:
    type: File[]
    doc: "Expression pathway analysis results"
    outputSource: differential_expression/pathway_results
  
  expression_pathway_plots:
    type: Directory
    doc: "Visualizations of differential expression analysis"
    outputSource: differential_expression/pathway_plots
  
  # 4. Methylation Analysis Results
  methylation_results:
    type: File[]
    doc: "Differential methylation analysis results"
    outputSource: methylation_analysis/methylation_results
  
  methylation_pathway_results:
    type: File
    doc: "Methylation pathway analysis results"
    outputSource: methylation_analysis/pathway_results
  
  methylation_pathway_plots:
    type: Directory
    doc: "Visualizations of methylation analysis"
    outputSource: methylation_analysis/pathway_plots
  
  # 5. Post-Clustering Analysis Results
  comparison_results:
    type: File
    doc: "Statistical comparison results between clusters"
    outputSource: post_clustering/comparison_results
    
  heatmaps:
    type: Directory
    doc: "Heatmap visualizations"
    outputSource: post_clustering/heatmaps
    
  bubble_plots:
    type: Directory
    doc: "Bubble plot visualizations"
    outputSource: post_clustering/bubble_plots
    
  survival_plots:
    type: Directory
    doc: "Survival analysis plots"
    outputSource: post_clustering/survival_plots
    
  sankey_plots:
    type: Directory
    doc: "Sankey diagram visualizations"
    outputSource: post_clustering/sankey_plots