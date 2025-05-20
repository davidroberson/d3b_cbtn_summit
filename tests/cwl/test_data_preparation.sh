#!/bin/bash

# Test script for the data_preparation tool
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
echo -e "${GREEN}      Testing Data Preparation Tool                               ${NC}"
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
TOOL_PATH="${WORKFLOW_DIR}/tools/data_preparation.cwl"
TEST_DATA_DIR="${SCRIPT_DIR}/test_data"
TEST_INPUTS="${TEST_DATA_DIR}/test_inputs.yaml"
OUTPUT_DIR="${SCRIPT_DIR}/test_output/data_preparation"

# Create output directory if it doesn't exist
echo -e "${YELLOW}Creating output directory: ${OUTPUT_DIR}${NC}"
mkdir -p "${OUTPUT_DIR}"

# Update the test input paths to match the current directory structure
echo -e "${YELLOW}Creating temporary test inputs file...${NC}"
TMP_INPUTS="${OUTPUT_DIR}/tmp_inputs.yaml"
cat > "${TMP_INPUTS}" << EOF
histology_file:
  class: File
  path: ${TEST_DATA_DIR}/histologies.tsv
short_histology: "HGAT"
count_file:
  class: File
  path: ${TEST_DATA_DIR}/gene-counts-rsem-expected_count-collapsed.rds
methyl_file:
  class: File
  path: ${TEST_DATA_DIR}/methyl-beta-values.rds
splice_file:
  class: File
  path: ${TEST_DATA_DIR}/psi-se-primary.func.rds
gtf_file:
  class: File
  path: ${TEST_DATA_DIR}/gencode.v39.primary_assembly.annotation.gtf.gz
num_features: 100
EOF

# Verify test data exists
echo -e "${YELLOW}Verifying test data files...${NC}"
for file in $(grep -o "${TEST_DATA_DIR}.*\.(tsv|rds|gz)" "${TMP_INPUTS}" | tr -d '"'); do
    if [ -f "$file" ]; then
        echo -e "  ${GREEN}✓${NC} $file"
    else
        echo -e "  ${RED}✗${NC} $file (not found!)"
        EXIT_CODE=1
    fi
done

if [ -n "$EXIT_CODE" ]; then
    echo -e "${RED}Error: Some test data files are missing!${NC}"
    exit 1
fi

# Print information about the test
echo -e "${YELLOW}Tool: ${TOOL_PATH}${NC}"
echo -e "${YELLOW}Test inputs: ${TMP_INPUTS}${NC}"
echo

# Run the tool
echo -e "${GREEN}Starting data preparation tool...${NC}"
echo

cwltool --force-docker-pull --outdir "${OUTPUT_DIR}" "${TOOL_PATH}" "${TMP_INPUTS}"

# Check if tool completed successfully
if [ $? -eq 0 ]; then
    echo
    echo -e "${GREEN}==================================================================${NC}"
    echo -e "${GREEN}      Data preparation tool completed successfully!              ${NC}"
    echo -e "${GREEN}      Results are available in: ${OUTPUT_DIR}                    ${NC}"
    echo -e "${GREEN}==================================================================${NC}"
else
    echo
    echo -e "${RED}==================================================================${NC}"
    echo -e "${RED}      Data preparation tool execution failed!                     ${NC}"
    echo -e "${RED}==================================================================${NC}"
    exit 1
fi

# Cleanup
echo -e "${YELLOW}Cleaning up temporary files...${NC}"
rm -f "${TMP_INPUTS}"