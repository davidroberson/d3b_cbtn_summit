# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!require("reshape2", quietly = TRUE))
  install.packages("reshape2")

if (!require("msigdbr", quietly = TRUE))
  install.packages("msigdbr")

if (!require("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr")

# Bioc packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!require("rtracklayer", quietly = TRUE))
  BiocManager::install("rtracklayer")

if (!require("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
