# Package Analysis for Bioconductor Workflow

## Packages Already in Bioconductor Base Image
- survival
- ggplot2 (via tidyverse)
- optparse
- gridExtra
- reshape2

## Packages We Need to Add
1. **Bioconductor packages:**
   - DESeq2
   - limma
   - clusterProfiler
   - rtracklayer
   - ComplexHeatmap
   - missMethyl
   - IlluminaHumanMethylationEPICanno.ilm10b4.hg19

2. **CRAN packages:**
   - IntNMF
   - datawizard
   - ggpubr
   - corrplot
   - circlize
   - mclust
   - msigdbr
   - survminer
   - tidyverse (even though parts might exist)

## Minimal Docker Image to Create
For an optimized approach, we could create a Dockerfile using Bioconductor as base:

```dockerfile
FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install required CRAN packages
RUN R -e 'install.packages(c("IntNMF", "datawizard", "ggpubr", "corrplot", "circlize", "mclust", "msigdbr", "survminer", "tidyverse"), repos="https://cran.rstudio.com/")'

# Install required Bioconductor packages
RUN R -e 'BiocManager::install(c("DESeq2", "limma", "clusterProfiler", "rtracklayer", "ComplexHeatmap", "missMethyl", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"))'

# Create working directory
WORKDIR /data
```

This would significantly improve build time compared to starting from scratch, while ensuring all required packages are available.