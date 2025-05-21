FROM pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.2

# Install additional R packages
RUN R -e "install.packages('IntNMF', dependencies = TRUE)"
RUN R -e "install.packages('survminer', repos='https://cran.rstudio.com/', dependencies = TRUE)"

# Create directory structure
RUN mkdir -p /app/cbtn_multiomics

# Copy R scripts and utilities from src directory
COPY src/cbtn_multiomics /app/cbtn_multiomics

# Set working directory
WORKDIR /app
