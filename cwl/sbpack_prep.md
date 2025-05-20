# CWL Package Preparation for sbpack

This document explains the changes made to the CWL tools to make them compatible with sbpack, which is used to upload CWL workflows to Seven Bridges platforms.

## Changes to CWL Tools

We've embedded the R scripts directly into the CWL files instead of using `$include` directives. This approach offers several advantages:

1. **Standalone files**: Each CWL file is self-contained with all required code
2. **Simpler packaging**: No need to resolve external references during packaging
3. **Guaranteed compatibility**: Works with all CWL tools without special handling

## Example: Data Preparation Tool

The `data_preparation.cwl` tool embeds the R script directly:

```yaml
baseCommand: ["Rscript", "--vanilla"]

requirements:
  - class: DockerRequirement
    dockerPull: "bioconductor/bioconductor_docker:RELEASE_3_18"
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.R
        entry: |
          #!/usr/bin/env Rscript
          # prepare files for multi-modal clustering
          suppressPackageStartupMessages({
            library(optparse)
            library(tidyverse)
            library(datawizard)
            library(reshape2)
            library(rtracklayer)
          })
          
          # ... rest of the script ...
          
      - entryname: utils/filter_cnv.R
        entry: |
          # filter cnvs by oncogenes/tsgs as well as copy number status cutoff
          filter_cnv <- function(myCNVData, myCancerGenes = cancer_genes) {
            # ... function implementation ...
          }
```

## Example: IntNMF Tool

Similar approach for the `integrative_nmf.cwl` tool:

```yaml
baseCommand: ["Rscript", "--vanilla"]

requirements:
  - class: DockerRequirement
    dockerPull: "bioconductor/bioconductor_docker:RELEASE_3_18"
  - class: InitialWorkDirRequirement
    listing:
      - entryname: script.R
        entry: |
          #!/usr/bin/env Rscript
          # run multi-modal clustering
          suppressPackageStartupMessages({
            library(optparse)
            library(IntNMF)
            library(tidyverse)
          })
          
          # ... rest of the script ...
          
      - entryname: utils/run_clusterstats.R
        entry: |
          # This script is to calculate cluster stats for each k
          suppressPackageStartupMessages({
            library(IntNMF)
            library(tidyverse)
          })
          
          run_clusterstats <- function(dat, output_dir, k_value) {
            # ... function implementation ...
          }
```

## Using sbpack

To package and upload the workflow to a Seven Bridges platform:

1. Install sbpack:
   ```
   pip install sbpack
   ```

2. Package and upload the workflow:
   ```
   sbpack platform project_id workflow.cwl
   ```

sbpack will automatically:
- Package all the CWL files with their embedded scripts
- Upload everything to the specified project

This approach ensures that all code dependencies are properly included in the uploaded workflow, making it fully portable and reproducible on the Seven Bridges platform.