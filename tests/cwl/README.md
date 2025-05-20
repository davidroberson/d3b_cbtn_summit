# CWL Workflow Tests

This directory contains test scripts for the CBTN Multi-Omic Clustering workflow and its individual tools. The tests ensure that each component of the workflow functions correctly with the provided test data.

## Test Data

Test data files are located in the `test_data/` directory. These are small, synthetic datasets suitable for testing the workflow functionality without requiring large computational resources.

## Test Scripts

### Individual Tool Tests

Each CWL tool has a dedicated test script:

- **`test_data_preparation.sh`**: Tests the data preparation step, which filters and transforms input data
- **`test_integrative_nmf.sh`**: Tests the integrative clustering step using the IntNMF algorithm
- **`test_dge_analysis.sh`**: Tests the differential gene expression analysis
- **`test_methylation_analysis.sh`**: Tests the methylation analysis step
- **`test_post_clustering.sh`**: Tests the post-clustering associations and visualizations

### Complete Workflow Test

- **`test_workflow.sh`**: Tests the entire workflow end-to-end

## Running Tests

### Prerequisites

- `cwltool` must be installed (`pip install cwltool`)
- Docker must be installed and running
- Bioconductor Docker image will be pulled automatically

### Running Individual Tool Tests

Each tool can be tested independently:

```bash
./test_data_preparation.sh
./test_integrative_nmf.sh
./test_dge_analysis.sh
./test_methylation_analysis.sh
./test_post_clustering.sh
```

Note that some tests depend on outputs from previous steps. For example, the integrative_nmf test requires outputs from the data_preparation test. The scripts will automatically run prerequisites if needed.

### Running the Complete Workflow

To test the entire workflow at once:

```bash
./test_workflow.sh
```

This will run all steps of the workflow in sequence.

## Test Outputs

Test outputs are saved in the `test_output/` directory, with a subdirectory for each tool:

- `test_output/data_preparation/`: Data preparation outputs
- `test_output/integrative_nmf/`: Clustering outputs
- `test_output/dge_analysis/`: Differential gene expression analysis outputs
- `test_output/methylation_analysis/`: Methylation analysis outputs
- `test_output/post_clustering/`: Post-clustering analysis outputs
- `test_output/workflow/`: Complete workflow outputs

## Troubleshooting

If tests fail, check the following:

1. Ensure Docker is running
2. Verify that cwltool is installed correctly
3. Check that test data files exist and are accessible
4. Examine the error messages for specific issues
5. Make sure you have adequate disk space and memory