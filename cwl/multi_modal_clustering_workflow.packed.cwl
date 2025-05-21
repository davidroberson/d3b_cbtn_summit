class: Workflow
cwlVersion: v1.2
doc: 'A workflow for multi-modal clustering analysis of cancer genomics data.

  Takes RNA-seq expression, methylation, and alternative splicing data

  as input to generate integrated clusters and perform downstream analyses.

  '
inputs:
- {doc: Clinical and histology data for samples (.tsv), id: histology_file, type: File}
- {doc: 'Short histology code to filter samples (e.g., ''HGAT'')', id: short_histology,
  type: string}
- {doc: RNA-seq expression matrix (.rds), id: count_file, type: File}
- {doc: Methylation beta values matrix (.rds), id: methyl_file, type: File}
- {doc: Alternative splice junction PSI values (.rds), id: splice_file, type: File}
- {doc: Gene annotation GTF file (.gtf.gz), id: gtf_file, type: File}
- {default: intNMF, doc: Clustering method to use (currently only intNMF supported),
  id: cluster_method, type: string}
- {default: 1000, doc: Number of top variable features to select for each data modality,
  id: num_features, type: int}
- {default: 10, doc: Maximum number of clusters to try, id: max_k, type: int}
label: CBTN Multi-Omic Clustering Workflow
outputs:
- {doc: Processed RNA-seq data matrix, id: processed_rna_data, outputSource: data_preparation/rna_data,
  type: File}
- {doc: Processed methylation data matrix, id: processed_methyl_data, outputSource: data_preparation/methyl_data,
  type: File}
- {doc: Processed splicing data matrix, id: processed_splice_data, outputSource: data_preparation/splice_data,
  type: File}
- {doc: Sample ID mapping file, id: sample_mapping, outputSource: data_preparation/samples_map,
  type: File}
- {doc: Cluster assignments for each sample, id: cluster_assignments, outputSource: cluster_analysis/cluster_assignments,
  type: File}
- {doc: RNA feature importance scores, id: feature_scores_rna, outputSource: cluster_analysis/feature_scores_rna,
  type: File}
- {doc: Methylation feature importance scores, id: feature_scores_methyl, outputSource: cluster_analysis/feature_scores_methyl,
  type: File}
- {doc: Splicing feature importance scores, id: feature_scores_splice, outputSource: cluster_analysis/feature_scores_splice,
  type: File}
- {doc: Consensus matrix plot for clustering, id: consensus_plot, outputSource: cluster_analysis/consensus_plot,
  type: File}
- {doc: Silhouette plot for clustering, id: silhouette_plot, outputSource: cluster_analysis/silhouette_plot,
  type: File}
- {doc: Differential expression analysis results, id: deseq_results, outputSource: differential_expression/deseq_results,
  type: File}
- doc: Expression pathway analysis results
  id: expression_pathway_results
  outputSource: differential_expression/pathway_results
  pickValue: all_non_null
  type: {items: File, type: array}
- {doc: Visualizations of differential expression analysis, id: expression_pathway_plots,
  outputSource: differential_expression/pathway_plots, type: Directory}
- doc: Differential methylation analysis results
  id: methylation_results
  outputSource: methylation_analysis/methylation_results
  pickValue: all_non_null
  type: {items: File, type: array}
- {doc: Methylation pathway analysis results, id: methylation_pathway_results, outputSource: methylation_analysis/pathway_results,
  type: File}
- {doc: Visualizations of methylation analysis, id: methylation_pathway_plots, outputSource: methylation_analysis/pathway_plots,
  type: Directory}
- {doc: Statistical comparison results between clusters, id: comparison_results, outputSource: post_clustering/comparison_results,
  type: File}
- {doc: Heatmap visualizations, id: heatmaps, outputSource: post_clustering/heatmaps,
  type: Directory}
- {doc: Bubble plot visualizations, id: bubble_plots, outputSource: post_clustering/bubble_plots,
  type: Directory}
- {doc: Survival analysis plots, id: survival_plots, outputSource: post_clustering/survival_plots,
  type: Directory}
- {doc: Sankey diagram visualizations, id: sankey_plots, outputSource: post_clustering/sankey_plots,
  type: Directory}
requirements:
- {class: SubworkflowFeatureRequirement}
- {class: ScatterFeatureRequirement}
- {class: MultipleInputFeatureRequirement}
- {class: InlineJavascriptRequirement}
steps:
- id: data_preparation
  in:
  - {id: histology_file, source: histology_file}
  - {id: short_histology, source: short_histology}
  - {id: count_file, source: count_file}
  - {id: methyl_file, source: methyl_file}
  - {id: splice_file, source: splice_file}
  - {id: gtf_file, source: gtf_file}
  - {id: num_features, source: num_features}
  out: [rna_data, methyl_data, splice_data, samples_map]
  run:
    baseCommand: [Rscript, --vanilla, /app/cbtn_multiomics/data_preparation/01-multi-modal-clustering-prepare-data.R]
    class: CommandLineTool
    cwlVersion: v1.2
    doc: 'Prepares RNA-seq, methylation, and alternative splicing data

      for multi-modal clustering analysis. Filters and transforms data

      according to specified parameters.

      '
    inputs:
    - id: histology_file
      inputBinding: {position: 1, prefix: --histology_file}
      type: File
    - default: HGAT
      id: short_histology
      inputBinding: {position: 2, prefix: --short_histology}
      type: string
    - id: count_file
      inputBinding: {position: 3, prefix: --count_file}
      type: File
    - id: methyl_file
      inputBinding: {position: 4, prefix: --methyl_file}
      type: File
    - id: splice_file
      inputBinding: {position: 5, prefix: --splice_file}
      type: File
    - id: gtf_file
      inputBinding: {position: 6, prefix: --gtf_file}
      type: File
    - default: 1000
      id: num_features
      inputBinding: {position: 7, prefix: --num_features}
      type: int
    - default: results
      id: output_dir
      inputBinding: {position: 8, prefix: --output_dir}
      type: string
    label: Data Preparation for Multi-Omic Clustering
    outputs:
    - id: rna_data
      outputBinding: {glob: $(inputs.output_dir)/rna_data.tsv}
      type: File
    - id: methyl_data
      outputBinding: {glob: $(inputs.output_dir)/methyl_data.tsv}
      type: File
    - id: splice_data
      outputBinding: {glob: $(inputs.output_dir)/splice_data.tsv}
      type: File
    - id: samples_map
      outputBinding: {glob: $(inputs.output_dir)/samples_map.tsv}
      type: File
    requirements:
    - {class: DockerRequirement, dockerPull: 'pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3'}
    - {class: InlineJavascriptRequirement}
    sbg:appVersion: [v1.2]
- id: cluster_analysis
  in:
  - {id: samples_map, source: data_preparation/samples_map}
  - {id: rna_data, source: data_preparation/rna_data}
  - {id: methyl_data, source: data_preparation/methyl_data}
  - {id: splice_data, source: data_preparation/splice_data}
  - {id: max_k, source: max_k}
  - {id: method, source: cluster_method}
  out: [cluster_assignments, feature_scores_rna, feature_scores_methyl, feature_scores_splice,
    consensus_plot, silhouette_plot]
  run:
    baseCommand: [Rscript, --vanilla, /app/cbtn_multiomics/integrative_nmf/01-multi-modal-clustering-run.R]
    class: CommandLineTool
    cwlVersion: v1.2
    doc: 'Performs multi-modal clustering using the Integrative Non-negative Matrix

      Factorization (IntNMF) method on multiple data types.

      '
    inputs:
    - id: samples_map
      inputBinding: {position: 1, prefix: --samples_map}
      type: File
    - id: rna_data
      inputBinding: {position: 2, prefix: --rna_file}
      type: File
    - id: methyl_data
      inputBinding: {position: 3, prefix: --methyl_file}
      type: File
    - id: splice_data
      inputBinding: {position: 4, prefix: --splice_file}
      type: File
    - default: 10
      id: max_k
      inputBinding: {position: 5, prefix: --max_k}
      type: int
    - default: intNMF
      id: method
      inputBinding: {position: 6, prefix: --cluster_method}
      type: string
    - default: results
      id: results_dir
      inputBinding: {position: 7, prefix: --results_dir}
      type: string
    - default: plots
      id: plots_dir
      inputBinding: {position: 8, prefix: --plots_dir}
      type: string
    label: Integrative NMF Clustering
    outputs:
    - id: cluster_assignments
      outputBinding: {glob: $(inputs.results_dir)/intnmf_clusters.tsv}
      type: File
    - id: feature_scores_rna
      outputBinding: {glob: $(inputs.results_dir)/feature_scores/feature_scores_rna.tsv}
      type: File
    - id: feature_scores_methyl
      outputBinding: {glob: $(inputs.results_dir)/feature_scores/feature_scores_methyl.tsv}
      type: File
    - id: feature_scores_splice
      outputBinding: {glob: $(inputs.results_dir)/feature_scores/feature_scores_splice.tsv}
      type: File
    - id: consensus_plot
      outputBinding: {glob: $(inputs.plots_dir)/intnmf_consensus_plot.pdf}
      type: File
    - id: silhouette_plot
      outputBinding: {glob: $(inputs.plots_dir)/intnmf_silhouette_plot.pdf}
      type: File
    requirements:
    - {class: DockerRequirement, dockerPull: 'pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3'}
    - {class: InlineJavascriptRequirement}
    sbg:appVersion: [v1.2]
- id: differential_expression
  in:
  - {id: expr_mat, source: count_file}
  - {id: gtf_file, source: gtf_file}
  - {id: cluster_file, source: cluster_analysis/cluster_assignments}
  out: [deseq_results, pathway_results, pathway_plots]
  run:
    baseCommand: [Rscript, --vanilla, /app/cbtn_multiomics/dge_analysis/01-dge_analysis_deseq.R]
    class: CommandLineTool
    cwlVersion: v1.2
    doc: 'Performs differential gene expression analysis between clusters

      and the rest of the samples using DESeq2, and conducts pathway

      enrichment analysis on the results.

      '
    inputs:
    - id: expr_mat
      inputBinding: {position: 1, prefix: --expr_mat}
      type: File
    - id: gtf_file
      inputBinding: {position: 2, prefix: --gtf_file}
      type: File
    - id: cluster_file
      inputBinding: {position: 3, prefix: --cluster_file}
      type: File
    - default: results
      id: results_dir
      inputBinding: {position: 4, prefix: --results_dir}
      type: string
    - default: plots
      id: plots_dir
      inputBinding: {position: 5, prefix: --plots_dir}
      type: string
    label: Differential Gene Expression Analysis
    outputs:
    - id: deseq_results
      outputBinding: {glob: $(inputs.results_dir)/diffexpr_output_per_cluster.tsv}
      type: File
    - id: pathway_results
      outputBinding:
        glob: [$(inputs.results_dir)/hallmark/*.tsv, $(inputs.results_dir)/reactome/*.tsv,
          '*.tsv_*']
      type: {items: File, type: array}
    - id: pathway_plots
      outputBinding: {glob: $(inputs.plots_dir)}
      type: Directory
    requirements:
    - {class: DockerRequirement, dockerPull: 'pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3'}
    - {class: InlineJavascriptRequirement}
    sbg:appVersion: [v1.2]
- id: methylation_analysis
  in:
  - {id: methyl_file, source: methyl_file}
  - {id: cluster_file, source: cluster_analysis/cluster_assignments}
  out: [methylation_results, pathway_results, pathway_plots]
  run:
    baseCommand: [Rscript, --vanilla, /app/cbtn_multiomics/methylation/01-limma_analysis.R]
    class: CommandLineTool
    cwlVersion: v1.2
    doc: 'Performs differential methylation analysis between clusters

      and conducts pathway enrichment analysis on the differentially

      methylated sites.

      '
    inputs:
    - id: methyl_file
      inputBinding: {position: 1, prefix: --methyl_file}
      type: File
    - id: cluster_file
      inputBinding: {position: 2, prefix: --cluster_file}
      type: File
    - default: results
      id: results_dir
      inputBinding: {position: 3, prefix: --results_dir}
      type: string
    - default: plots
      id: plots_dir
      inputBinding: {position: 4, prefix: --plots_dir}
      type: string
    label: Methylation Analysis
    outputs:
    - id: methylation_results
      outputBinding: {glob: $(inputs.results_dir)/limma_output/*.tsv}
      type: {items: File, type: array}
    - id: pathway_results
      outputBinding: {glob: $(inputs.results_dir)/dms_gsameth_output/hallmark/genebody_promoter_gsameth_output_per_cluster.tsv}
      type: File
    - id: pathway_plots
      outputBinding: {glob: $(inputs.plots_dir)}
      type: Directory
    requirements:
    - {class: DockerRequirement, dockerPull: 'pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3'}
    - {class: InlineJavascriptRequirement}
    sbg:appVersion: [v1.2]
- id: post_clustering
  in:
  - {id: cluster_file, source: cluster_analysis/cluster_assignments}
  - {id: histology_file, source: histology_file}
  - {id: rna_data, source: data_preparation/rna_data}
  - {id: methyl_data, source: data_preparation/methyl_data}
  - {id: splice_data, source: data_preparation/splice_data}
  - {id: feature_scores_rna, source: cluster_analysis/feature_scores_rna}
  - {id: feature_scores_methyl, source: cluster_analysis/feature_scores_methyl}
  - {id: feature_scores_splice, source: cluster_analysis/feature_scores_splice}
  out: [comparison_results, heatmaps, bubble_plots, survival_plots, sankey_plots]
  run:
    arguments:
    - {position: 0, valueFrom: "set -e\n\n# Set environment variables for run_analysis.sh\nexport
        cluster_file=$(inputs.cluster_file.path)\nexport histology_file=$(inputs.histology_file.path)\nexport
        count_file=$(inputs.rna_data.path)\nexport methyl_file=$(inputs.methyl_data.path)\nexport
        splice_file=$(inputs.splice_data.path)\nexport feature_scores_rna=$(inputs.feature_scores_rna.path)\nexport
        feature_scores_methyl=$(inputs.feature_scores_methyl.path)\nexport feature_scores_splice=$(inputs.feature_scores_splice.path)\n\n#
        Create output directories\nmkdir -p results/intnmf/survival_os\nmkdir -p results/intnmf/survival_efs\nmkdir
        -p plots/intnmf/heatmaps\nmkdir -p plots/intnmf/bubble_plots\nmkdir -p plots/intnmf/survival_os\nmkdir
        -p plots/intnmf/survival_efs\nmkdir -p plots/intnmf/sankey_plots\n\n# Create
        the expected folder structure for the adapted R scripts\nmkdir -p plots/heatmaps\nmkdir
        -p plots/bubble_plots\nmkdir -p plots/survival_plots\nmkdir -p plots/sankey_plots\n\n#
        Run the analysis script\nbash /app/cbtn_multiomics/post_clustering/run_analysis.sh\n\n#
        Copy results to the expected output paths for CWL\ncp -r plots/intnmf/heatmaps/*
        plots/heatmaps/ || true\ncp -r plots/intnmf/bubble_plots/* plots/bubble_plots/
        || true\ncp -r plots/intnmf/survival_os/* plots/survival_plots/ || true\ncp
        -r plots/intnmf/survival_efs/* plots/survival_plots/ || true\ncp -r plots/intnmf/sankey_plots/*
        plots/sankey_plots/ || true\n\n# Make sure we have the correct output file
        for comparison_results\nif [ -f \"results/intnmf/ari_mm_vs_subtypes.txt\"
        ] && [ -f \"results/intnmf/chisq_mm_vs_subtypes.txt\" ]; then\n  # Combine
        them into the file expected by the CWL\n  cat \"results/intnmf/ari_mm_vs_subtypes.txt\"
        > \"results/intnmf/chisq_ari_mm_vs_subtypes.txt\"\n  echo -e \"\\n\\n\" >>
        \"results/intnmf/chisq_ari_mm_vs_subtypes.txt\"\n  cat \"results/intnmf/chisq_mm_vs_subtypes.txt\"
        >> \"results/intnmf/chisq_ari_mm_vs_subtypes.txt\"\nfi\n\n# Remove placeholder
        files\nfind plots -name \"placeholder.txt\" -delete\n"}
    baseCommand: [/bin/bash, -c]
    class: CommandLineTool
    cwlVersion: v1.2
    doc: 'Performs post-clustering analysis including statistical comparisons,

      visualization, and survival analysis.

      '
    inputs:
    - {id: cluster_file, type: File}
    - {id: histology_file, type: File}
    - {id: rna_data, type: File}
    - {id: methyl_data, type: File}
    - {id: splice_data, type: File}
    - {id: feature_scores_rna, type: File}
    - {id: feature_scores_methyl, type: File}
    - {id: feature_scores_splice, type: File}
    label: Post-Clustering Analysis
    outputs:
    - id: comparison_results
      outputBinding: {glob: results/intnmf/chisq_ari_mm_vs_subtypes.txt}
      type: File
    - id: heatmaps
      outputBinding: {glob: plots/heatmaps}
      type: Directory
    - id: bubble_plots
      outputBinding: {glob: plots/bubble_plots}
      type: Directory
    - id: survival_plots
      outputBinding: {glob: plots/survival_plots}
      type: Directory
    - id: sankey_plots
      outputBinding: {glob: plots/sankey_plots}
      type: Directory
    requirements:
    - {class: DockerRequirement, dockerPull: 'pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3'}
    - {class: InlineJavascriptRequirement}
    - class: EnvVarRequirement
      envDef: {data_dir: /tmp, output_dir: results/intnmf, plots_dir: plots/intnmf}
    sbg:appVersion: [v1.2]
