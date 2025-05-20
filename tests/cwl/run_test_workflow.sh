#!/bin/bash

# Script to run the multi-omic clustering workflow with test data
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
echo -e "${GREEN}      Running Multi-Omic Clustering Workflow with Test Data      ${NC}"
echo -e "${GREEN}==================================================================${NC}"

# Check if cwltool is installed
if ! command -v cwltool &> /dev/null; then
    echo -e "${RED}Error: cwltool is not installed${NC}"
    echo "Please install it with: pip install cwltool"
    exit 1
fi

# Define paths
WORKFLOW_DIR=$(dirname "$(readlink -f "$0")")
WORKFLOW_PATH="${WORKFLOW_DIR}/multi_modal_clustering_workflow.cwl"
TEST_INPUTS="${WORKFLOW_DIR}/test_data/test_inputs.yaml"
OUTPUT_DIR="${WORKFLOW_DIR}/test_output"

# Create output directory if it doesn't exist
echo -e "${YELLOW}Creating output directory: ${OUTPUT_DIR}${NC}"
mkdir -p "${OUTPUT_DIR}"

# Verify test data exists
echo -e "${YELLOW}Verifying test data files...${NC}"
for file in $(grep -o "/workspaces.*\.(tsv|rds|gz)" "${TEST_INPUTS}" | tr -d '"'); do
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

# Print information about the workflow
echo -e "${YELLOW}Workflow file: ${WORKFLOW_PATH}${NC}"
echo -e "${YELLOW}Test inputs: ${TEST_INPUTS}${NC}"
echo

# Run the workflow
echo -e "${GREEN}Starting workflow execution...${NC}"
echo

cwltool --outdir "${OUTPUT_DIR}" "${WORKFLOW_PATH}" "${TEST_INPUTS}"

# Check if workflow completed successfully
if [ $? -eq 0 ]; then
    echo
    echo -e "${GREEN}==================================================================${NC}"
    echo -e "${GREEN}      Workflow completed successfully!                          ${NC}"
    echo -e "${GREEN}      Results are available in: ${OUTPUT_DIR}                   ${NC}"
    echo -e "${GREEN}==================================================================${NC}"
else
    echo
    echo -e "${RED}==================================================================${NC}"
    echo -e "${RED}      Workflow execution failed!                                ${NC}"
    echo -e "${RED}==================================================================${NC}"
    exit 1
fi