#!/usr/bin/env cwl-runner

cwlVersion: v1.2
class: CommandLineTool

label: Data Preparation for Multi-Omic Clustering
doc: |
  Prepares RNA-seq, methylation, and alternative splicing data
  for multi-modal clustering analysis. Filters and transforms data
  according to specified parameters.

baseCommand: ["Rscript", "--vanilla", "/app/cbtn_multiomics/data_preparation/01-multi-modal-clustering-prepare-data.R"]

requirements:
  - class: DockerRequirement
    # Use our custom image with all required packages and scripts pre-installed
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.3"
  - class: InlineJavascriptRequirement

inputs:
  histology_file:
    type: File
    inputBinding:
      position: 1
      prefix: --histology_file
  
  short_histology:
    type: string
    default: "HGAT"
    inputBinding:
      position: 2
      prefix: --short_histology
  
  count_file:
    type: File
    inputBinding:
      position: 3
      prefix: --count_file
  
  methyl_file:
    type: File
    inputBinding:
      position: 4
      prefix: --methyl_file
  
  splice_file:
    type: File
    inputBinding:
      position: 5
      prefix: --splice_file
  
  gtf_file:
    type: File
    inputBinding:
      position: 6
      prefix: --gtf_file
  
  num_features:
    type: int
    default: 1000
    inputBinding:
      position: 7
      prefix: --num_features
  
  output_dir:
    type: string
    default: "results"
    inputBinding:
      position: 8
      prefix: --output_dir

outputs:
  rna_data:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/rna_data.tsv
  
  methyl_data:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/methyl_data.tsv
  
  splice_data:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/splice_data.tsv
  
  samples_map:
    type: File
    outputBinding:
      glob: $(inputs.output_dir)/samples_map.tsv