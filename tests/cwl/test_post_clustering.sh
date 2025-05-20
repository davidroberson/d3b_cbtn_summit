#!/bin/bash

# Test script for the post_clustering tool
# Author: Claude
# Date: May 20, 2025

set -e  # Exit on error
set -o pipefail

# Define colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Print header
echo -e "${GREEN}==================================================================${NC}"
echo -e "${GREEN}      Testing Post-Clustering Analysis Tool                       ${NC}"
echo -e "${GREEN}==================================================================${NC}"

# Check if cwltool is installed
if ! command -v cwltool &> /dev/null; then
    echo -e "${RED}Error: cwltool is not installed${NC}"
    echo "Please install it with: pip install cwltool"
    exit 1
fi

# Define paths
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
WORKFLOW_DIR=$(realpath "${SCRIPT_DIR}/../../cwl")
TOOL_PATH="${WORKFLOW_DIR}/tools/post_clustering.cwl"
TEST_DATA_DIR="${SCRIPT_DIR}/test_data"
DATA_PREP_OUTPUT="${SCRIPT_DIR}/test_output/data_preparation"
CLUSTER_OUTPUT="${SCRIPT_DIR}/test_output/integrative_nmf"
OUTPUT_DIR="${SCRIPT_DIR}/test_output/post_clustering"

# Check if previous step outputs exist, run prerequisite tests if needed
if [ ! -f "${CLUSTER_OUTPUT}/cluster_assignments" ] || [ ! -f "${CLUSTER_OUTPUT}/feature_scores_rna" ] || [ ! -f "${CLUSTER_OUTPUT}/feature_scores_methyl" ] || [ ! -f "${CLUSTER_OUTPUT}/feature_scores_splice" ]; then
    echo -e "${YELLOW}Clustering outputs not found. Running integrative NMF test first...${NC}"
    "${SCRIPT_DIR}/test_integrative_nmf.sh"
    if [ $? -ne 0 ]; then
        echo -e "${RED}Integrative NMF test failed, cannot proceed with post-clustering analysis test.${NC}"
        exit 1
    fi
fi

if [ ! -f "${DATA_PREP_OUTPUT}/rna_data" ] || [ ! -f "${DATA_PREP_OUTPUT}/methyl_data" ] || [ ! -f "${DATA_PREP_OUTPUT}/splice_data" ]; then
    echo -e "${YELLOW}Data preparation outputs not found. Running data preparation test first...${NC}"
    "${SCRIPT_DIR}/test_data_preparation.sh"
    if [ $? -ne 0 ]; then
        echo -e "${RED}Data preparation test failed, cannot proceed with post-clustering analysis test.${NC}"
        exit 1
    fi
fi

# Create output directory if it doesn't exist
echo -e "${YELLOW}Creating output directory: ${OUTPUT_DIR}${NC}"
mkdir -p "${OUTPUT_DIR}"

# Create temporary inputs file for the post_clustering tool
echo -e "${YELLOW}Creating temporary test inputs file...${NC}"
TMP_INPUTS="${OUTPUT_DIR}/tmp_inputs.yaml"
cat > "${TMP_INPUTS}" << EOF
cluster_file:
  class: File
  path: ${CLUSTER_OUTPUT}/cluster_assignments
histology_file:
  class: File
  path: ${TEST_DATA_DIR}/histologies.tsv
rna_data:
  class: File
  path: ${DATA_PREP_OUTPUT}/rna_data
methyl_data:
  class: File
  path: ${DATA_PREP_OUTPUT}/methyl_data
splice_data:
  class: File
  path: ${DATA_PREP_OUTPUT}/splice_data
feature_scores_rna:
  class: File
  path: ${CLUSTER_OUTPUT}/feature_scores_rna
feature_scores_methyl:
  class: File
  path: ${CLUSTER_OUTPUT}/feature_scores_methyl
feature_scores_splice:
  class: File
  path: ${CLUSTER_OUTPUT}/feature_scores_splice
EOF

# Verify input files exist
echo -e "${YELLOW}Verifying input files...${NC}"
for file in $(grep -o "\(${TEST_DATA_DIR}\|${DATA_PREP_OUTPUT}\|${CLUSTER_OUTPUT}\)/.*" "${TMP_INPUTS}" | tr -d '"'); do
    if [ -f "$file" ]; then
        echo -e "  ${GREEN}✓${NC} $file"
    else
        echo -e "  ${RED}✗${NC} $file (not found!)"
        EXIT_CODE=1
    fi
done

if [ -n "$EXIT_CODE" ]; then
    echo -e "${RED}Error: Some input files are missing!${NC}"
    exit 1
fi

# Print information about the test
echo -e "${YELLOW}Tool: ${TOOL_PATH}${NC}"
echo -e "${YELLOW}Test inputs: ${TMP_INPUTS}${NC}"
echo

# Run the tool
echo -e "${GREEN}Starting post-clustering analysis tool...${NC}"
echo

cwltool --outdir "${OUTPUT_DIR}" "${TOOL_PATH}" "${TMP_INPUTS}"

# Check if tool completed successfully
if [ $? -eq 0 ]; then
    echo
    echo -e "${GREEN}==================================================================${NC}"
    echo -e "${GREEN}      Post-clustering analysis tool completed successfully!      ${NC}"
    echo -e "${GREEN}      Results are available in: ${OUTPUT_DIR}                    ${NC}"
    echo -e "${GREEN}==================================================================${NC}"
else
    echo
    echo -e "${RED}==================================================================${NC}"
    echo -e "${RED}      Post-clustering analysis tool execution failed!            ${NC}"
    echo -e "${RED}==================================================================${NC}"
    exit 1
fi

# Cleanup
echo -e "${YELLOW}Cleaning up temporary files...${NC}"
rm -f "${TMP_INPUTS}"