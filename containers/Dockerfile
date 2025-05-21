FROM bioconductor/bioconductor_docker:RELEASE_3_18

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install basic R packages
RUN R -e "install.packages('optparse', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('tidyverse', repos='https://cran.rstudio.com/', dependencies = TRUE)"

# Data preparation packages
RUN R -e "install.packages('datawizard', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('reshape2', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('matrixStats', repos='https://cran.rstudio.com/', dependencies = TRUE)"

# Clustering packages
RUN R -e "install.packages('cluster', repos='https://cran.rstudio.com/', dependencies = TRUE)"

# Visualization packages
#RUN R -e "install.packages('ggplot2', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('ggpubr', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('corrplot', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('gridExtra', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('circlize', repos='https://cran.rstudio.com/', dependencies = TRUE)"

# Statistical analysis packages
RUN R -e "install.packages('survival', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('mclust', repos='https://cran.rstudio.com/', dependencies = TRUE)"
RUN R -e "install.packages('msigdbr', repos='https://cran.rstudio.com/', dependencies = TRUE)"

# Bioconductor packages - core
RUN R -e "BiocManager::install('rtracklayer', ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('DESeq2', ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('limma', ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('ComplexHeatmap', ask = FALSE, update = FALSE)"

# Bioconductor packages - genomics
RUN R -e "BiocManager::install('clusterProfiler', ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('missMethyl', ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('IlluminaHumanMethylationEPICanno.ilm10b4.hg19', ask = FALSE, update = FALSE)"
