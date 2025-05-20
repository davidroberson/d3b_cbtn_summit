#!/bin/bash

# Combined test script for the multi-omic clustering workflow and its components
# Author: Claude
# Date: May 20, 2025

set -e  # Exit on error
set -o pipefail

# Define colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Define paths
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
WORKFLOW_DIR=$(realpath "${SCRIPT_DIR}/../../cwl")
TEST_DATA_DIR="${SCRIPT_DIR}/test_data"
OUTPUT_BASE_DIR="${SCRIPT_DIR}/test_output"

# Define step names for mapping and printing
declare -A step_names
step_names[1]="Data Preparation"
step_names[2]="Integrative NMF Clustering"
step_names[3]="Differential Gene Expression Analysis"
step_names[4]="Methylation Analysis"
step_names[5]="Post-Clustering Analysis"

# Define step tool paths
declare -A step_tools
step_tools[1]="${WORKFLOW_DIR}/tools/data_preparation.cwl"
step_tools[2]="${WORKFLOW_DIR}/tools/integrative_nmf.cwl"
step_tools[3]="${WORKFLOW_DIR}/tools/dge_analysis.cwl"
step_tools[4]="${WORKFLOW_DIR}/tools/methylation_analysis.cwl"
step_tools[5]="${WORKFLOW_DIR}/tools/post_clustering.cwl"

# Define step output directories
declare -A step_output_dirs
step_output_dirs[1]="${OUTPUT_BASE_DIR}/data_preparation"
step_output_dirs[2]="${OUTPUT_BASE_DIR}/integrative_nmf"
step_output_dirs[3]="${OUTPUT_BASE_DIR}/dge_analysis"
step_output_dirs[4]="${OUTPUT_BASE_DIR}/methylation_analysis"
step_output_dirs[5]="${OUTPUT_BASE_DIR}/post_clustering"

# Print usage information
function show_usage() {
    echo -e "${BOLD}Usage:${NC} $0 [OPTIONS]"
    echo
    echo -e "${BOLD}Description:${NC}"
    echo "  Run tests for the multi-omic clustering workflow components."
    echo
    echo -e "${BOLD}Options:${NC}"
    echo "  -h, --help         Show this help message"
    echo "  -a, --all          Run all tests in sequence (default)"
    echo "  -s, --step NUMBER  Run a specific step:"
    echo "                     1 = Data Preparation"
    echo "                     2 = Integrative NMF Clustering"
    echo "                     3 = Differential Gene Expression Analysis"
    echo "                     4 = Methylation Analysis"
    echo "                     5 = Post-Clustering Analysis"
    echo "  -c, --clean        Clean output directories before running tests"
    echo "  -v, --validate     Only validate the CWL files (no execution)"
    echo
    echo -e "${BOLD}Examples:${NC}"
    echo "  $0 --all           # Run all steps in sequence"
    echo "  $0 --step 2        # Run only the Integrative NMF step"
    echo "  $0 --step 1,3,5    # Run steps 1, 3, and 5"
    echo "  $0 --clean --all   # Clean output directories and run all steps"
}

# Validate CWL files
function validate_cwl() {
    echo -e "${YELLOW}Validating CWL files...${NC}"
    
    for i in {1..5}; do
        local tool="${step_tools[$i]}"
        echo -e "Validating ${step_names[$i]} (${tool})..."
        
        cwltool --validate "${tool}"
        if [ $? -eq 0 ]; then
            echo -e "  ${GREEN}✓${NC} ${step_names[$i]} validation passed"
        else
            echo -e "  ${RED}✗${NC} ${step_names[$i]} validation failed"
            EXIT_CODE=1
        fi
    done
    
    # Also validate the main workflow
    echo -e "Validating main workflow (${WORKFLOW_DIR}/multi_modal_clustering_workflow.cwl)..."
    cwltool --validate "${WORKFLOW_DIR}/multi_modal_clustering_workflow.cwl"
    if [ $? -eq 0 ]; then
        echo -e "  ${GREEN}✓${NC} Main workflow validation passed"
    else
        echo -e "  ${RED}✗${NC} Main workflow validation failed"
        EXIT_CODE=1
    fi
    
    if [ -n "$EXIT_CODE" ]; then
        echo -e "${RED}Validation failed for some CWL files.${NC}"
        exit 1
    else
        echo -e "${GREEN}All CWL files validated successfully.${NC}"
    fi
}

# Check prerequisites
function check_prerequisites() {
    # Check if cwltool is installed
    if ! command -v cwltool &> /dev/null; then
        echo -e "${RED}Error: cwltool is not installed${NC}"
        echo "Please install it with: pip install cwltool"
        exit 1
    fi
}

# Clean output directories
function clean_outputs() {
    echo -e "${YELLOW}Cleaning output directories...${NC}"
    
    for i in {1..5}; do
        local dir="${step_output_dirs[$i]}"
        if [ -d "$dir" ]; then
            echo -e "Removing ${dir}..."
            rm -rf "$dir"
        fi
    done
    
    echo -e "${GREEN}Output directories cleaned.${NC}"
}

# Verify test data exists
function verify_test_data() {
    echo -e "${YELLOW}Verifying test data files...${NC}"
    
    # List of required test data files
    local files=(
        "${TEST_DATA_DIR}/histologies.tsv"
        "${TEST_DATA_DIR}/gene-counts-rsem-expected_count-collapsed.rds"
        "${TEST_DATA_DIR}/methyl-beta-values.rds"
        "${TEST_DATA_DIR}/psi-se-primary.func.rds"
        "${TEST_DATA_DIR}/gencode.v39.primary_assembly.annotation.gtf.gz"
    )
    
    for file in "${files[@]}"; do
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
}

# Run a specific step
function run_step() {
    local step=$1
    local step_name="${step_names[$step]}"
    local tool="${step_tools[$step]}"
    local output_dir="${step_output_dirs[$step]}"
    
    echo -e "${GREEN}==================================================================${NC}"
    echo -e "${GREEN}      Testing ${step_name} (Step $step)                           ${NC}"
    echo -e "${GREEN}==================================================================${NC}"
    
    # Create output directory
    echo -e "${YELLOW}Creating output directory: ${output_dir}${NC}"
    mkdir -p "${output_dir}"
    
    # Create inputs for the specific step
    local tmp_inputs="${output_dir}/tmp_inputs.yaml"
    
    case $step in
        1) # Data Preparation
            cat > "${tmp_inputs}" << EOF
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
            ;;
            
        2) # Integrative NMF
            # Check if data preparation outputs exist
            if [ ! -f "${step_output_dirs[1]}/rna_data.tsv" ] || [ ! -f "${step_output_dirs[1]}/methyl_data.tsv" ] || [ ! -f "${step_output_dirs[1]}/splice_data.tsv" ] || [ ! -f "${step_output_dirs[1]}/samples_map.tsv" ]; then
                echo -e "${YELLOW}Data preparation outputs not found. Running data preparation first...${NC}"
                run_step 1
            fi
            
            cat > "${tmp_inputs}" << EOF
samples_map:
  class: File
  path: ${step_output_dirs[1]}/samples_map.tsv
rna_data:
  class: File
  path: ${step_output_dirs[1]}/rna_data.tsv
methyl_data:
  class: File
  path: ${step_output_dirs[1]}/methyl_data.tsv
splice_data:
  class: File
  path: ${step_output_dirs[1]}/splice_data.tsv
max_k: 3
method: "intNMF"
EOF
            ;;
            
        3) # DGE Analysis
            # Check if clustering outputs exist
            if [ ! -f "${step_output_dirs[2]}/intnmf_clusters.tsv" ]; then
                echo -e "${YELLOW}Clustering outputs not found. Running integrative NMF first...${NC}"
                run_step 2
            fi
            
            cat > "${tmp_inputs}" << EOF
expr_mat:
  class: File
  path: ${TEST_DATA_DIR}/gene-counts-rsem-expected_count-collapsed.rds
gtf_file:
  class: File
  path: ${TEST_DATA_DIR}/gencode.v39.primary_assembly.annotation.gtf.gz
cluster_file:
  class: File
  path: ${step_output_dirs[2]}/intnmf_clusters.tsv
histology_file:
  class: File
  path: ${TEST_DATA_DIR}/histologies.tsv
EOF
            ;;
            
        4) # Methylation Analysis
            # Check if clustering outputs exist
            if [ ! -f "${step_output_dirs[2]}/intnmf_clusters.tsv" ]; then
                echo -e "${YELLOW}Clustering outputs not found. Running integrative NMF first...${NC}"
                run_step 2
            fi
            
            cat > "${tmp_inputs}" << EOF
methyl_file:
  class: File
  path: ${TEST_DATA_DIR}/methyl-beta-values.rds
cluster_file:
  class: File
  path: ${step_output_dirs[2]}/intnmf_clusters.tsv
histology_file:
  class: File
  path: ${TEST_DATA_DIR}/histologies.tsv
EOF
            ;;
            
        5) # Post-Clustering Analysis
            # Check if previous step outputs exist
            if [ ! -f "${step_output_dirs[2]}/intnmf_clusters.tsv" ] || [ ! -f "${step_output_dirs[2]}/feature_scores_rna.tsv" ] || [ ! -f "${step_output_dirs[2]}/feature_scores_methyl.tsv" ] || [ ! -f "${step_output_dirs[2]}/feature_scores_splice.tsv" ]; then
                echo -e "${YELLOW}Clustering outputs not found. Running integrative NMF first...${NC}"
                run_step 2
            fi
            
            if [ ! -f "${step_output_dirs[1]}/rna_data.tsv" ] || [ ! -f "${step_output_dirs[1]}/methyl_data.tsv" ] || [ ! -f "${step_output_dirs[1]}/splice_data.tsv" ]; then
                echo -e "${YELLOW}Data preparation outputs not found. Running data preparation first...${NC}"
                run_step 1
            fi
            
            cat > "${tmp_inputs}" << EOF
cluster_file:
  class: File
  path: ${step_output_dirs[2]}/intnmf_clusters.tsv
histology_file:
  class: File
  path: ${TEST_DATA_DIR}/histologies.tsv
rna_data:
  class: File
  path: ${step_output_dirs[1]}/rna_data.tsv
methyl_data:
  class: File
  path: ${step_output_dirs[1]}/methyl_data.tsv
splice_data:
  class: File
  path: ${step_output_dirs[1]}/splice_data.tsv
feature_scores_rna:
  class: File
  path: ${step_output_dirs[2]}/feature_scores_rna.tsv
feature_scores_methyl:
  class: File
  path: ${step_output_dirs[2]}/feature_scores_methyl.tsv
feature_scores_splice:
  class: File
  path: ${step_output_dirs[2]}/feature_scores_splice.tsv
EOF
            ;;
    esac
    
    # Verify input files
    echo -e "${YELLOW}Verifying input files for ${step_name}...${NC}"
    for file in $(grep -o "path:.*" "${tmp_inputs}" | sed 's/path: //' | tr -d '"'); do
        if [ -f "$file" ]; then
            echo -e "  ${GREEN}✓${NC} $file"
        else
            echo -e "  ${RED}✗${NC} $file (not found!)"
            local input_error=1
        fi
    done
    
    if [ -n "$input_error" ]; then
        echo -e "${RED}Error: Some input files are missing!${NC}"
        exit 1
    fi
    
    # Print information about the test
    echo -e "${YELLOW}Tool: ${tool}${NC}"
    echo -e "${YELLOW}Test inputs: ${tmp_inputs}${NC}"
    echo
    
    # Run the tool
    echo -e "${GREEN}Starting ${step_name} tool...${NC}"
    echo
    
    cwltool --outdir "${output_dir}" "${tool}" "${tmp_inputs}"
    
    # Check if tool completed successfully
    if [ $? -eq 0 ]; then
        echo
        echo -e "${GREEN}==================================================================${NC}"
        echo -e "${GREEN}      ${step_name} tool completed successfully!                  ${NC}"
        echo -e "${GREEN}      Results are available in: ${output_dir}                    ${NC}"
        echo -e "${GREEN}==================================================================${NC}"
        
        # Clean up
        echo -e "${YELLOW}Cleaning up temporary files...${NC}"
        rm -f "${tmp_inputs}"
        
        return 0
    else
        echo
        echo -e "${RED}==================================================================${NC}"
        echo -e "${RED}      ${step_name} tool execution failed!                        ${NC}"
        echo -e "${RED}==================================================================${NC}"
        exit 1
    fi
}

# Run all steps in sequence
function run_all_steps() {
    echo -e "${GREEN}==================================================================${NC}"
    echo -e "${GREEN}      Running All Steps in Sequence                                ${NC}"
    echo -e "${GREEN}==================================================================${NC}"
    
    for i in {1..5}; do
        run_step $i
    done
    
    echo -e "${GREEN}==================================================================${NC}"
    echo -e "${GREEN}      All steps completed successfully!                           ${NC}"
    echo -e "${GREEN}==================================================================${NC}"
}

# Parse command line arguments
CLEAN=false
VALIDATE=false
RUN_ALL=true
STEPS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            show_usage
            exit 0
            ;;
        -a|--all)
            RUN_ALL=true
            shift
            ;;
        -s|--step)
            RUN_ALL=false
            STEPS="$2"
            shift 2
            ;;
        -c|--clean)
            CLEAN=true
            shift
            ;;
        -v|--validate)
            VALIDATE=true
            shift
            ;;
        *)
            echo -e "${RED}Error: Unknown option $1${NC}"
            show_usage
            exit 1
            ;;
    esac
done

# Check prerequisites
check_prerequisites

# Verify test data
verify_test_data

# Clean output if requested
if [ "$CLEAN" = true ]; then
    clean_outputs
fi

# Validate CWL files if requested
if [ "$VALIDATE" = true ]; then
    validate_cwl
    exit 0
fi

# Run the requested steps
if [ "$RUN_ALL" = true ]; then
    run_all_steps
else
    # Run specific steps
    IFS=',' read -ra STEP_ARRAY <<< "$STEPS"
    for step in "${STEP_ARRAY[@]}"; do
        if [[ "$step" =~ ^[1-5]$ ]]; then
            run_step "$step"
        else
            echo -e "${RED}Error: Invalid step number: $step (must be 1-5)${NC}"
            show_usage
            exit 1
        fi
    done
fi

echo -e "${GREEN}Done!${NC}"
exit 0