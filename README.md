# D3b CBTN Summit Workshop

## CBTN Multi-Omic Subtyping

This repository provides the necessary setup and analysis scripts to demonstrate the use of CBTN data for multi-omic clustering within the Cavatica Cloud Environment.

## CWL Workflow Implementation

We've implemented the multi-omic clustering workflow as a portable, reproducible workflow using the Common Workflow Language (CWL). The complete workflow is available in the [cwl directory](cwl/).

### Key Features

- **Portable**: Run the workflow locally, on HPC, or in the cloud
- **Reproducible**: Docker containers ensure consistent execution environments
- **Modular**: Individual steps can be reused in other workflows
- **Self-contained**: Includes test data and comprehensive documentation

### Getting Started with the CWL Workflow

```bash
# Clone the repository
git clone https://github.com/d3b-center/d3b-cbtn-summit.git
cd d3b-cbtn-summit/cwl

# Generate test data
make generate-test-data

# Run the workflow with test data
make test
```

See the [cwl/README.md](cwl/README.md) file for complete documentation of the CWL implementation.

## Original R Scripts

The original R scripts used in the workshop are available in the `analyses` directory:

1. `data_preparation/`: Prepares multiomics data for clustering
2. `intNMF/`: Runs integrative NMF clustering
3. `post_clustering_associations/`: Generates visualizations and statistical analyses
4. `dge_pathway_analysis/`: Performs differential gene expression and pathway analysis
5. `methylation_analysis/`: Performs methylation analysis

### Workshop Steps

In this workshop, we perform the following:
1. Open a new cloud-based RStudio instance within Cavatica's Data Studio
2. Clone the d3b-cbtn-summit github repository using the terminal within RStudio
3. Review the modules within the repository and their key functionality. 
4. Prepare epigenomic, transcriptomic, and alternatively splice transcriptome data for high-grade astrocytomas (HGATs) as input to multi-omic clustering. 
5. Run intNMF matrix factorization for multi-omic clustering using the three data layers outlined in (4). 
6. Evaluate novel subtypes in relation to known disease subtypes and survival characteristics. 
7. Evaluate cluster-specific differential expression patterns.
8. Time permitting, evaluate cluster-specific differential methylation and differentially methylated pathways.