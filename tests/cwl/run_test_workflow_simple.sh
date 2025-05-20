#!/bin/bash

# Simple script to run the multi-omic clustering workflow with test data using cwltool
# Author: Claude
# Date: May 20, 2025

set -e  # Exit on error
set -o pipefail

# Define workflow and data paths
WORKFLOW_DIR=$(dirname "$(readlink -f "$0")")
WORKFLOW_PATH="${WORKFLOW_DIR}/multi_modal_clustering_workflow.cwl"
TEST_DATA_DIR="${WORKFLOW_DIR}/test_data"
OUTPUT_DIR="${WORKFLOW_DIR}/test_output"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "Running workflow with cwltool and Docker..."
echo

# Run cwltool with direct command line arguments instead of YAML
cwltool \
  --outdir "${OUTPUT_DIR}" \
  --enable-dev \
  --no-read-only \
  --preserve-environment PATH HOME \
  --docker \
  "${WORKFLOW_PATH}" \
  --histology_file "${TEST_DATA_DIR}/histologies.tsv" \
  --short_histology "HGAT" \
  --count_file "${TEST_DATA_DIR}/gene-counts-rsem-expected_count-collapsed.rds" \
  --methyl_file "${TEST_DATA_DIR}/methyl-beta-values.rds" \
  --splice_file "${TEST_DATA_DIR}/psi-se-primary.func.rds" \
  --gtf_file "${TEST_DATA_DIR}/gencode.v39.primary_assembly.annotation.gtf.gz" \
  --num_features 100 \
  --max_k 3

exit_code=$?
if [ $exit_code -eq 0 ]; then
    echo
    echo "Workflow completed successfully!"
    echo "Results are available in: ${OUTPUT_DIR}"
else
    echo
    echo "Workflow execution failed with exit code: ${exit_code}"
    exit 1
fi