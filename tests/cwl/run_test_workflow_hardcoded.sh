#!/bin/bash

# Super simple script with hardcoded paths to run the multi-omic clustering workflow
# Author: Claude
# Date: May 20, 2025

# Create output directory
mkdir -p /workspaces/d3b_cbtn_summit/cwl/test_output

# Check if cwltool is installed
if ! command -v cwltool &> /dev/null; then
    echo "Error: cwltool is not installed"
    echo "Please install it with: pip install cwltool"
    exit 1
fi

# Print the command for clarity
echo "Running the following command:"
echo "cwltool \\"
echo "  --outdir /workspaces/d3b_cbtn_summit/cwl/test_output \\"
echo "  --enable-dev \\"
echo "  --docker \\"
echo "  /workspaces/d3b_cbtn_summit/cwl/multi_modal_clustering_workflow.cwl \\"
echo "  --histology_file /workspaces/d3b_cbtn_summit/cwl/test_data/histologies.tsv \\"
echo "  --short_histology HGAT \\"
echo "  --count_file /workspaces/d3b_cbtn_summit/cwl/test_data/gene-counts-rsem-expected_count-collapsed.rds \\"
echo "  --methyl_file /workspaces/d3b_cbtn_summit/cwl/test_data/methyl-beta-values.rds \\"
echo "  --splice_file /workspaces/d3b_cbtn_summit/cwl/test_data/psi-se-primary.func.rds \\"
echo "  --gtf_file /workspaces/d3b_cbtn_summit/cwl/test_data/gencode.v39.primary_assembly.annotation.gtf.gz \\"
echo "  --num_features 100 \\"
echo "  --max_k 3"
echo

# Run cwltool with direct command line arguments and hardcoded paths
cwltool \
  --outdir /workspaces/d3b_cbtn_summit/cwl/test_output \
  --enable-dev \
  --force-docker-pull \
  /workspaces/d3b_cbtn_summit/cwl/multi_modal_clustering_workflow.cwl \
  --histology_file /workspaces/d3b_cbtn_summit/cwl/test_data/histologies.tsv \
  --short_histology HGAT \
  --count_file /workspaces/d3b_cbtn_summit/cwl/test_data/gene-counts-rsem-expected_count-collapsed.rds \
  --methyl_file /workspaces/d3b_cbtn_summit/cwl/test_data/methyl-beta-values.rds \
  --splice_file /workspaces/d3b_cbtn_summit/cwl/test_data/psi-se-primary.func.rds \
  --gtf_file /workspaces/d3b_cbtn_summit/cwl/test_data/gencode.v39.primary_assembly.annotation.gtf.gz \
  --num_features 100 \
  --max_k 3