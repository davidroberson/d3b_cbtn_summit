# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!require("datawizard", quietly = TRUE))
  install.packages("datawizard")

if (!require("reshape2", quietly = TRUE))
  install.packages("reshape2")

# Bioc packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("rtracklayer", quietly = TRUE))
  BiocManager::install("rtracklayer")
