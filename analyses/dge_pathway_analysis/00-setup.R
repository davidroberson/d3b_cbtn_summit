# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", repos = "http://cran.us.r-project.org")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")

if (!require("reshape2", quietly = TRUE))
  install.packages("reshape2", repos = "http://cran.us.r-project.org")

if (!require("msigdbr", quietly = TRUE))
  install.packages("msigdbr", repos = "http://cran.us.r-project.org")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2", repos = "http://cran.us.r-project.org")

if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr", repos = "http://cran.us.r-project.org")

# Bioc packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")

if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!require("rtracklayer", quietly = TRUE))
  BiocManager::install("rtracklayer")

if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
