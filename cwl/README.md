# Building and Pushing the Docker Image

This directory contains files needed to build the Docker image for the multi-modal clustering workflow.

## Docker Image Requirements

The workflow requires a Docker image with the following R packages:
- Core bioinformatics: DESeq2, limma, clusterProfiler, missMethyl
- Clustering: IntNMF
- Utilities: optparse, tidyverse, reshape2
- Visualization: ComplexHeatmap, ggplot2, survminer

## Building the Docker Image

Two Dockerfile options are provided:

1. **Dockerfile** - A full build from r-base (time-intensive)
2. **Dockerfile.simple** - A simpler build extending the Bioconductor Docker image (recommended)

### Recommended Approach

Use the `build_push_docker.sh` script to build and push the Docker image to the CAVATICA registry:

```bash
# Make the script executable if needed
chmod +x build_push_docker.sh

# Run the script
./build_push_docker.sh
```

This will:
1. Build the Docker image using Dockerfile.simple
2. Tag it as pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0
3. Push it to the CAVATICA Docker registry

### Manual Build

If you prefer to build manually:

```bash
# Build the image
docker build -f Dockerfile.simple -t pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0 .

# Push to CAVATICA
docker push pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0
```

## Updating the CWL Files

The CWL tool files already reference this Docker image. If you change the image name or tag, you'll need to update the DockerRequirement in each tool file:

```yaml
requirements:
  - class: DockerRequirement
    dockerPull: "pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.0"
```

## Troubleshooting

If you encounter issues:

1. **Authentication errors**: Ensure you're logged in to the CAVATICA Docker registry:
   ```
   docker login pgc-images.sbgenomics.com
   ```

2. **Install errors**: Try installing R packages in smaller groups or one at a time.

3. **Package compatibility issues**: Check for version conflicts and adjust the Dockerfile accordingly.# Multi-Omic Clustering Workflow

This directory contains a Common Workflow Language (CWL) implementation of a multi-omic clustering pipeline for cancer genomics data.

## Implementation Approach

We have chosen an approach that uses the official Bioconductor Docker image with dynamic package installation:

1. **Base Container**: All tools use the `bioconductor/bioconductor_docker:RELEASE_3_18` image
2. **Dynamic Package Installation**: Each R script includes code to install required packages if they're not already available
3. **No Custom Docker Needed**: This eliminates the need to build and maintain a custom Docker image

### Benefits of this approach:

- **Simplified Deployment**: No need to build and push custom Docker images
- **Reproducibility**: Uses a stable, well-maintained Bioconductor image
- **Flexibility**: Easy to add or update R packages without rebuilding containers
- **Maintenance**: Less overhead in managing Docker images

## Running the Workflow

### Input Parameters

Prepare your input parameters file (example in `inputs.yaml`):

```yaml
histology_file:
  class: File
  path: /path/to/histologies.tsv
short_histology: "HGAT"  # Histology type to filter
count_file:
  class: File
  path: /path/to/gene-counts-rsem-expected_count-collapsed.rds
methyl_file:
  class: File
  path: /path/to/methyl-beta-values.rds
splice_file:
  class: File
  path: /path/to/psi-se-primary.func.rds
gtf_file:
  class: File
  path: /path/to/gencode.v39.primary_assembly.annotation.gtf.gz
num_features: 1000  # Number of top variable features to select
max_k: 10  # Maximum number of clusters to try
```

### Running on CAVATICA

The workflow has been uploaded to CAVATICA and can be run using the platform's web interface or API.

Project: `david.roberson/copy-of-cbtn-summit-2024`
Workflow: `multi-modal-clustering-bioc`

## Workflow Structure

The workflow consists of five main steps:

1. **Data Preparation**: Filters and transforms input data matrices
2. **Integrative NMF Clustering**: Performs multi-modal clustering
3. **Differential Expression Analysis**: Identifies gene expression patterns in clusters
4. **Methylation Analysis**: Identifies differential methylation patterns
5. **Post-Clustering Analysis**: Generates visualizations and statistics

## Performance Considerations

- **First Run Time**: The first run may take longer as packages are installed dynamically
- **Subsequent Runs**: If using the same execution environment, packages will remain installed
- **Memory Requirements**: Some steps (especially clustering) require significant memory

## Troubleshooting

- **Package Installation Failures**: Check error messages for missing system dependencies
- **Memory Errors**: Try increasing available memory for the task
- **File Format Issues**: Ensure input files match expected formats (RDS files, TSV, etc.)

## Future Improvements

- Add test datasets
- Add parameter validation
- Implement resource optimization# Running the Multi-Omic Clustering Workflow on CAVATICA

The workflow has been successfully uploaded to CAVATICA and is ready to use\!

## Current Implementation

Our workflow uses the official Bioconductor Docker image () which has many common packages pre-installed. The R scripts have been modified to dynamically install any missing packages at runtime.

## Requirements for Running on CAVATICA

When running on CAVATICA, ensure that:

1. Instance type has sufficient memory (at least 16GB recommended)
2. Container has write permissions to install packages
3. Input files are properly uploaded to CAVATICA storage

## Package Analysis

After analyzing the Bioconductor base image, we found these are the additional packages needed:

### Bioconductor packages:
- DESeq2
- limma
- clusterProfiler
- rtracklayer
- ComplexHeatmap
- missMethyl
- IlluminaHumanMethylationEPICanno.ilm10b4.hg19

### CRAN packages:
- IntNMF
- datawizard
- ggpubr
- corrplot
- circlize
- mclust
- msigdbr
- survminer
- tidyverse (parts)

## Test Data

We've created small test datasets that demonstrate the workflow functionality:
- Small matrices (1000 features x 20 samples)
- Simplified GTF annotation
- Generated random clinical data
- Reduced clustering parameters

## Future Enhancements

1. Build a custom Docker image with all packages pre-installed
2. Add more comprehensive test datasets
3. Implement workflow-level parameter validation
4. Add parallelization options for larger datasets

## Contact

For any issues or questions about this workflow, please contact your workflow administrator.
