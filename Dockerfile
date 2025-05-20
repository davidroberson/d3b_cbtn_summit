FROM pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering:v1.0.2


RUN R -e "install.packages('IntNMF', dependencies = TRUE)"
