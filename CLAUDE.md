# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository contains R scripts for multi-omic clustering analysis of CBTN (Children's Brain Tumor Network) data. It provides a complete workflow for analyzing and clustering high-grade astrocytomas (HGATs) using epigenomic, transcriptomic, and alternatively spliced transcriptome data.

The repository is structured with two parallel implementations:
1. R scripts in the `analyses/` directory for direct execution
2. A portable Common Workflow Language (CWL) implementation in the `cwl/` directory

## Analysis Workflow

The workflow consists of multiple analysis modules, executed in the following order:

1. **Data Preparation** (`/analyses/data_preparation/`): Filters and transforms RNA, methylation, and splicing data for HGAT samples.
2. **IntNMF Clustering** (`/analyses/intNMF/`): Performs multi-modal clustering using the IntNMF algorithm.
3. **Post-Clustering Associations** (`/analyses/post_clustering_associations/`): Generates survival plots, heatmaps, and association analyses.
4. **Differential Gene Expression Analysis** (`/analyses/dge_pathway_analysis/`): Performs differential gene expression analysis and pathway enrichment.
5. **Methylation Analysis** (`/analyses/methylation_analysis/`): Performs differential methylation analysis and pathway enrichment.

## Common Commands

### Running R Analysis Modules

#### Running a Complete Analysis Module

Each analysis module can be run using its respective `run_analysis.sh` script:

```bash
cd analyses/data_preparation
bash run_analysis.sh
```

Similar commands apply for other modules:

```bash
cd analyses/intNMF
bash run_analysis.sh
```

```bash
cd analyses/post_clustering_associations
bash run_analysis.sh
```

```bash
cd analyses/dge_pathway_analysis
bash run_analysis.sh
```

```bash
cd analyses/methylation_analysis
bash run_analysis.sh
```

#### Running Individual Analysis Steps

Each module's R scripts can be run individually with appropriate parameters. For example:

```bash
cd analyses/data_preparation
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
cd analyses/[module_name]
Rscript --vanilla 00-setup.R
```

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

#### Building and Pushing the Docker Image

The CWL workflow uses a Docker container with all required packages pre-installed:

```bash
# Build and push the Docker image to CAVATICA
./utils/docker/build_push_docker.sh

# Or build manually
docker build -f Dockerfile.simple -t pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0 .
docker push pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0
```

#### Uploading the Workflow to CAVATICA

To upload the workflow and test data to CAVATICA:

```bash
# Using the all-in-one Python script
python utils/upload_scripts/upload_to_cavatica_all_in_one.py
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

#### Bioconductor Packages
- DESeq2: Differential gene expression analysis
- rtracklayer: Genomic interval manipulation
- clusterProfiler: Pathway enrichment analysis

## Data Inputs

The analysis requires multiple input files:
- Histology data for HGAT samples
- RNA-seq expression data (expected counts)
- Methylation beta values
- Alternative splicing (PSI) data
- Gencode GTF annotation

The main data directory is set to `/sbgenomics/project-files/opc-v15/v15/` by default in the scripts.

## Expected Outputs

Each module produces its own outputs:
- Filtered and transformed data matrices
- Clustering results and feature scores
- Statistical comparisons and visualizations
- Differential expression/methylation results
- Pathway enrichment results

The outputs are saved in `results/` and `plots/` directories within each module.

## Project Structure

### R Analysis Modules (`analyses/`)

- `data_preparation/`: Filters and transforms RNA, methylation, and splicing data
- `intNMF/`: Performs multi-modal clustering using the IntNMF algorithm
- `post_clustering_associations/`: Generates visualizations and statistical analyses
- `dge_pathway_analysis/`: Performs differential gene expression and pathway analysis
- `methylation_analysis/`: Performs methylation analysis

### CWL Workflow (`cwl/`)

- `multi_modal_clustering_workflow.cwl`: Main workflow definition
- `tools/`: Individual tool definitions for each analysis step
- `test_data/`: Small test datasets and scripts to create them
- `docs/`: Workflow documentation and diagrams

### Utility Scripts (`utils/`)

- `docker/build_push_docker.sh`: Builds and pushes the Docker image
- `upload_scripts/`: Various Python scripts for uploading workflows to CAVATICA

### Tests (`tests/`)

- `cwl/`: Test scripts and data for CWL workflow
- `unit/`: Unit tests for individual components
- `integration/`: Integration tests for the full workflow

## Troubleshooting

### Common Issues with R Scripts

- **Memory Errors**: The default memory limit is set to 100GB in the `run_analysis.sh` scripts with `echo "R_MAX_VSIZE=100Gb" > .Renviron`
- **Missing Packages**: Each module's `00-setup.R` script installs required packages

### Common Issues with CWL

- **Docker Authentication**: Ensure you're logged in to the CAVATICA Docker registry with `docker login pgc-images.sbgenomics.com`
- **Package Installation Failures**: Check error messages for missing system dependencies
- **Memory Errors**: Try increasing available memory for resource-intensive steps