# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This repository contains R scripts for multi-omic clustering analysis of CBTN (Children's Brain Tumor Network) data. It provides a complete workflow for analyzing and clustering high-grade astrocytomas (HGATs) using epigenomic, transcriptomic, and alternatively spliced transcriptome data.

## Analysis Workflow

The workflow consists of multiple analysis modules, executed in the following order:

1. **Data Preparation** (`/analyses/data_preparation/`): Filters and transforms RNA, methylation, and splicing data for HGAT samples.
2. **IntNMF Clustering** (`/analyses/intNMF/`): Performs multi-modal clustering using the IntNMF algorithm.
3. **Post-Clustering Associations** (`/analyses/post_clustering_associations/`): Generates survival plots, heatmaps, and association analyses.
4. **Differential Gene Expression Analysis** (`/analyses/dge_pathway_analysis/`): Performs differential gene expression analysis and pathway enrichment.
5. **Methylation Analysis** (`/analyses/methylation_analysis/`): Performs differential methylation analysis and pathway enrichment.

## Common Commands

### Running a Complete Analysis Module

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

### Running Individual Analysis Steps

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

### Setting Up the R Environment

Each module has a setup script that installs the necessary R packages:

```bash
cd analyses/[module_name]
Rscript --vanilla 00-setup.R
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