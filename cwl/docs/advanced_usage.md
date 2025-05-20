# Advanced Usage Guide

This document provides advanced usage instructions for the multi-omic clustering workflow.

## Parameter Customization

### Clustering Parameters

The following parameters can be adjusted to optimize clustering results:

- **num_features**: Number of top variable features to select from each data modality
  - Default: 1000
  - Recommendation: 500-2000 depending on dataset size
  - Increasing this value improves feature coverage but increases memory usage

- **max_k**: Maximum number of clusters to try
  - Default: 10
  - Recommendation: Start with 3-12, then refine based on silhouette scores
  - Higher values increase computation time exponentially

### Performance Optimization

For large datasets, consider these strategies:

1. **Memory Management**:
   - Reduce `num_features` to lower memory requirements
   - Use cloud instances with sufficient memory (16GB+ recommended)
   - Monitor memory usage during execution

2. **Disk Space**:
   - Ensure at least 50GB free space for intermediate files
   - Clean up temporary files between runs

3. **Runtime Optimization**:
   - Split execution into stages for better error recovery
   - Use pre-filtered datasets when possible

## Custom Container Creation

For reproducible environments with all packages pre-installed:

1. **Using Provided Dockerfiles**:

```bash
# Full image with all dependencies
cd cwl/containers
docker build -t my-org/cbtn-multiomic:full -f Dockerfile .

# Simplified build using Bioconductor base
docker build -t my-org/cbtn-multiomic:simple -f Dockerfile.simple .
```

2. **Customizing Container**:

```dockerfile
# Example: Add custom packages to Dockerfile.simple
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install additional R packages
RUN R -e 'BiocManager::install(c("yourPackage1", "yourPackage2"))'

# Add custom scripts
COPY ./scripts /opt/scripts

# Set working directory
WORKDIR /data
```

3. **Updating CWL Files to Use Custom Container**:

Edit DockerRequirement sections in CWL tool definitions:

```yaml
requirements:
  - class: DockerRequirement
    dockerPull: "my-org/cbtn-multiomic:simple"
```

## Extended Analyses

### Integrating Additional Data Types

This workflow can be extended to integrate other data types:

1. **Copy Number Variation (CNV)**:
   - Create a new data preparation tool for CNV data
   - Add CNV data matrix to IntNMF input
   - Implement CNV-specific post-clustering analyses

2. **Proteomics**:
   - Follow similar patterns to RNA-seq data preparation
   - Adjust normalization methods as needed

### Custom Downstream Analyses

Additional downstream analyses can be added:

1. **Immune Cell Deconvolution**:
   ```R
   if (!require("CIBERSORT"))
     install.packages("CIBERSORT")
   
   # Add immune deconvolution analysis
   deconv_results <- CIBERSORT(expr_data, reference)
   ```

2. **Drug Response Prediction**:
   ```R
   if (!require("PharmacoGx"))
     BiocManager::install("PharmacoGx")
   
   # Predict drug response
   drug_response <- predictSensitivity(expr_data, model)
   ```

## Troubleshooting

### Common Issues

1. **Package Installation Failures**:
   - Error: `unable to install packages... permission denied`
   - Solution: Use `--user` flag when installing packages or ensure container has write access

2. **Memory Limitations**:
   - Error: `cannot allocate vector of size X Gb`
   - Solution: Increase available memory or reduce dataset size/feature count

3. **Missing Dependencies**:
   - Error: `there is no package called 'X'`
   - Solution: Check package installation section in scripts

### Runtime Debugging

For better visibility into workflow execution:

```bash
# Run with verbose logging
cwltool --debug multi_modal_clustering_workflow.cwl inputs.yaml

# Run individual steps for debugging
cwltool --no-container tools/data_preparation.cwl step_inputs.yaml
```

## Workflow Extension

To extend the workflow with new components:

1. **Create a new CWL tool definition**:
   ```yaml
   cwlVersion: v1.2
   class: CommandLineTool
   
   baseCommand: ["Rscript", "--vanilla"]
   
   requirements:
     - class: DockerRequirement
       dockerPull: "bioconductor/bioconductor_docker:RELEASE_3_18"
   
   inputs:
     new_input:
       type: File
       inputBinding:
         position: 1
         prefix: --input
   
   outputs:
     new_output:
       type: File
       outputBinding:
         glob: "results/output.tsv"
   ```

2. **Update the main workflow**:
   ```yaml
   steps:
     new_step:
       run: tools/new_tool.cwl
       in:
         new_input: previous_step/output
       out:
         - new_output
   ```

## CAVATICA-Specific Optimizations

For optimal performance on CAVATICA:

1. **Instance Selection**:
   - Use `c4.2xlarge` (8 vCPU/16GB) for most analyses
   - Consider `r4.2xlarge` (8 vCPU/32GB) for large datasets

2. **Storage Configuration**:
   - Mount temporary storage for large intermediate files
   - Use project storage for final outputs only

3. **App Configuration**:
   - Set appropriate instance hint flags
   - Configure retry parameters for network issues

4. **Batch Processing**:
   - Use CAVATICA batch functionality for cohort analyses
   - Structure batch inputs for parallel execution