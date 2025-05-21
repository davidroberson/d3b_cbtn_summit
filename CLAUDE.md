# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository contains R scripts for multi-omic clustering analysis of CBTN (Children's Brain Tumor Network) data. It provides a complete workflow for analyzing and clustering high-grade astrocytomas (HGATs) using epigenomic, transcriptomic, and alternatively spliced transcriptome data.

The repository is structured with two parallel implementations:
1. R scripts in the `src/cbtn_multiomics/` directory for direct execution
2. A portable Common Workflow Language (CWL) implementation in the `cwl/` directory

## Analysis Workflow

The workflow consists of multiple analysis modules, executed in the following order:

1. **Data Preparation** (`/src/cbtn_multiomics/data_preparation/`): Filters and transforms RNA, methylation, and splicing data for HGAT samples.
2. **IntNMF Clustering** (`/src/cbtn_multiomics/integrative_nmf/`): Performs multi-modal clustering using the IntNMF algorithm.
3. **Post-Clustering Associations** (`/src/cbtn_multiomics/post_clustering/`): Generates survival plots, heatmaps, and association analyses.
4. **Differential Gene Expression Analysis** (`/src/cbtn_multiomics/dge_analysis/`): Performs differential gene expression analysis and pathway enrichment.
5. **Methylation Analysis** (`/src/cbtn_multiomics/methylation/`): Performs differential methylation analysis and pathway enrichment.

## Common Commands

### Running R Analysis Modules

#### Running a Complete Analysis Module

Each analysis module can be run using its respective `run_analysis.sh` script:

```bash
cd src/cbtn_multiomics/data_preparation
bash run_analysis.sh
```

Similar commands apply for other modules:

```bash
cd src/cbtn_multiomics/integrative_nmf
bash run_analysis.sh
```

```bash
cd src/cbtn_multiomics/post_clustering
bash run_analysis.sh
```

```bash
cd src/cbtn_multiomics/dge_analysis
bash run_analysis.sh
```

```bash
cd src/cbtn_multiomics/methylation
bash run_analysis.sh
```

#### Running Individual Analysis Steps

Each module's R scripts can be run individually with appropriate parameters. For example:

```bash
cd src/cbtn_multiomics/data_preparation
Rscript --vanilla 01-multi-modal-clustering-prepare-data.R \
--histology_file $histology_file \
--short_histology "HGAT" \
--count_file $count_file \
--methyl_file $methyl_b_file \
--splice_file $functional_sites_splice_file \
--gtf_file $gtf_file \
--output_dir "results"
```

#### Setting Up the R Environment

Each module has a setup script that installs the necessary R packages:

```bash
cd src/cbtn_multiomics/[module_name]
Rscript --vanilla 00-setup.R
```

Note: The R scripts are designed to automatically install missing packages when they run, making them more robust across different environments. Each script includes checks for required packages and will install them if needed.

### Running the CWL Workflow

#### Running with Test Data

The CWL workflow comes with test data and a Makefile for easy execution:

```bash
cd cwl
make generate-test-data  # Creates small test datasets
make test                # Runs the workflow with test data
```

#### Running with Custom Data

To run the workflow with your own data:

```bash
cd cwl
# Edit inputs.yaml with your data paths
cwltool multi_modal_clustering_workflow.cwl inputs.yaml
```

#### Running Individual CWL Workflow Steps

You can run individual steps of the workflow for testing:

```bash
cd cwl
make test-data-preparation  # Test just the data preparation step
make test-clustering        # Test just the clustering step
```

#### Running CWL Tests

The repository contains test scripts for each component of the workflow:

```bash
cd tests/cwl
./run_all_tests.sh       # Run all tests
./test_data_preparation.sh    # Test only the data preparation step
./test_integrative_nmf.sh     # Test only the clustering step
./test_dge_analysis.sh        # Test only the DGE analysis step
./test_methylation_analysis.sh # Test only the methylation analysis step
./test_post_clustering.sh     # Test only the post-clustering analysis step
./test_workflow.sh           # Test the entire workflow end-to-end
```

#### Validating CWL Files

You can validate all CWL files in the workflow:

```bash
cd cwl
make validate
```

#### Building and Pushing the Docker Image

The CWL workflow uses a Docker container with all required packages pre-installed:

```bash
# Build and push the Docker image to CAVATICA
./utils/docker/build_push_docker.sh

# Or build manually
docker build -f containers/Dockerfile.simple -t pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0 .
docker push pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0
```

#### Uploading the Workflow to CAVATICA

To upload the workflow and test data to CAVATICA:

```bash
# Using the all-in-one Python script
python utils/upload_scripts/upload_to_cavatica_all_in_one.py

# Or using the Makefile
cd cwl
make upload-cavatica
```

#### Cleaning Up Temporary Files

After running tests or workflows, you can clean up temporary files:

```bash
cd cwl
make clean
```

#### Getting Help with the Makefile

To see all available commands in the Makefile:

```bash
cd cwl
make help
```

## Key Dependencies

### R Packages

#### CRAN Packages
- optparse: Command-line option parsing
- tidyverse: Data manipulation and visualization
- IntNMF: Integrative non-negative matrix factorization
- datawizard: Data transformations
- reshape2: Data reshaping
- ggplot2: Data visualization
- ggpubr: Publication-ready plots
- msigdbr: MSigDB gene sets
- corrplot: Correlation plots
- circlize: Circular visualization
- mclust: Model-based clustering
- survminer: Survival analysis visualization

#### Bioconductor Packages
- DESeq2: Differential gene expression analysis
- rtracklayer: Genomic interval manipulation
- clusterProfiler: Pathway enrichment analysis
- limma: Linear models for microarray data
- ComplexHeatmap: Complex heatmap visualization
- missMethyl: Methylation analysis
- IlluminaHumanMethylationEPICanno.ilm10b4.hg19: Methylation annotation

## Data Inputs

The analysis requires multiple input files:
- Histology data for HGAT samples
- RNA-seq expression data (expected counts)
- Methylation beta values
- Alternative splicing (PSI) data
- Gencode GTF annotation

The main data directory is set to `/sbgenomics/project-files/opc-v15/v15/` by default in the scripts, but paths can be modified in the run_analysis.sh files or provided as command-line arguments.

## Expected Outputs

Each module produces its own outputs:
- Filtered and transformed data matrices
- Clustering results and feature scores
- Statistical comparisons and visualizations
- Differential expression/methylation results
- Pathway enrichment results

The outputs are saved in `results/` and `plots/` directories within each module.

## Project Structure

### R Analysis Modules (`src/cbtn_multiomics/`)

- `data_preparation/`: Filters and transforms RNA, methylation, and splicing data
- `integrative_nmf/`: Performs multi-modal clustering using the IntNMF algorithm
- `post_clustering/`: Generates visualizations and statistical analyses
- `dge_analysis/`: Performs differential gene expression and pathway analysis
- `methylation/`: Performs methylation analysis

### CWL Workflow (`cwl/`)

- `multi_modal_clustering_workflow.cwl`: Main workflow definition
- `tools/`: Individual tool definitions for each analysis step
- `test_data/`: Small test datasets and scripts to create them
- `Makefile`: Commands for validating, testing, and deploying the workflow

### Tests (`tests/`)

- `cwl/`: Test scripts and data for CWL workflow
- `unit/`: Unit tests for individual components
- `integration/`: Integration tests for the full workflow

### Utility Scripts (`utils/`)

- `docker/build_push_docker.sh`: Builds and pushes the Docker image
- `upload_scripts/`: Various Python scripts for uploading workflows to CAVATICA

## Troubleshooting

### Common Issues with R Scripts

- **Memory Errors**: The default memory limit is set to 100GB in the `run_analysis.sh` scripts with `echo "R_MAX_VSIZE=100Gb" > .Renviron`
- **Missing Packages**: Each module's R scripts now include inline package installation checks and will automatically install missing packages
- **Data Format Issues**: Some scripts (especially methylation analysis) include error handling and fallback methods for problematic data formats

### Common Issues with CWL

- **Docker Authentication**: Ensure you're logged in to the CAVATICA Docker registry with `docker login pgc-images.sbgenomics.com`
- **Package Installation Failures**: Check error messages for missing system dependencies
- **Memory Errors**: Try increasing available memory for resource-intensive steps

## Recent Updates

- R scripts have been updated to include robust package installation checks and will automatically install missing dependencies
- Added error handling and data validation in analysis scripts, particularly in the methylation analysis module
- Added fallback functionality for handling various data format issues