#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Post-Clustering Analysis
doc: |
  Performs post-clustering analysis including statistical comparisons,
  visualization, and survival analysis.

baseCommand: ["/bin/bash", "-c"]

arguments:
  - position: 0
    valueFrom: |
      set -e
      
      # Set environment variables for run_analysis.sh
      export cluster_file=$(inputs.cluster_file.path)
      export histology_file=$(inputs.histology_file.path)
      export count_file=$(inputs.rna_data.path)
      export methyl_file=$(inputs.methyl_data.path)
      export splice_file=$(inputs.splice_data.path)
      export feature_scores_rna=$(inputs.feature_scores_rna.path)
      export feature_scores_methyl=$(inputs.feature_scores_methyl.path)
      export feature_scores_splice=$(inputs.feature_scores_splice.path)
      
      # Create output directories
      mkdir -p results/intnmf/survival_os
      mkdir -p results/intnmf/survival_efs
      mkdir -p plots/intnmf/heatmaps
      mkdir -p plots/intnmf/bubble_plots
      mkdir -p plots/intnmf/survival_os
      mkdir -p plots/intnmf/survival_efs
      mkdir -p plots/intnmf/sankey_plots
      
      # Create the expected folder structure for the adapted R scripts
      mkdir -p plots/heatmaps
      mkdir -p plots/bubble_plots
      mkdir -p plots/survival_plots
      mkdir -p plots/sankey_plots
      
      # Run the analysis script
      bash /app/cbtn_multiomics/post_clustering/run_analysis.sh
      
      # Copy results to the expected output paths for CWL
      cp -r plots/intnmf/heatmaps/* plots/heatmaps/ || true
      cp -r plots/intnmf/bubble_plots/* plots/bubble_plots/ || true
      cp -r plots/intnmf/survival_os/* plots/survival_plots/ || true
      cp -r plots/intnmf/survival_efs/* plots/survival_plots/ || true
      cp -r plots/intnmf/sankey_plots/* plots/sankey_plots/ || true
      
      # Make sure we have the correct output file for comparison_results
      if [ -f "results/intnmf/ari_mm_vs_subtypes.txt" ] && [ -f "results/intnmf/chisq_mm_vs_subtypes.txt" ]; then
        # Combine them into the file expected by the CWL
        cat "results/intnmf/ari_mm_vs_subtypes.txt" > "results/intnmf/chisq_ari_mm_vs_subtypes.txt"
        echo -e "\n\n" >> "results/intnmf/chisq_ari_mm_vs_subtypes.txt"
        cat "results/intnmf/chisq_mm_vs_subtypes.txt" >> "results/intnmf/chisq_ari_mm_vs_subtypes.txt"
      fi
      
      # Remove placeholder files
      find plots -name "placeholder.txt" -delete

requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3"
  - class: InlineJavascriptRequirement
  - class: EnvVarRequirement
    envDef:
      data_dir: "/tmp"
      output_dir: "results/intnmf"
      plots_dir: "plots/intnmf"

inputs:
  cluster_file:
    type: File
  
  histology_file:
    type: File
  
  rna_data:
    type: File
  
  methyl_data:
    type: File
  
  splice_data:
    type: File
  
  feature_scores_rna:
    type: File
  
  feature_scores_methyl:
    type: File
  
  feature_scores_splice:
    type: File

outputs:
  comparison_results:
    type: File
    outputBinding:
      glob: "results/intnmf/chisq_ari_mm_vs_subtypes.txt"
  
  heatmaps:
    type: Directory
    outputBinding:
      glob: "plots/heatmaps"
  
  bubble_plots:
    type: Directory
    outputBinding:
      glob: "plots/bubble_plots"
  
  survival_plots:
    type: Directory
    outputBinding:
      glob: "plots/survival_plots"
  
  sankey_plots:
    type: Directory
    outputBinding:
      glob: "plots/sankey_plots"