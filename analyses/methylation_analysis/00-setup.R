# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", repos = "http://cran.us.r-project.org")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")

if (!require("data.table", quietly = TRUE))
  install.packages("data.table", repos = "http://cran.us.r-project.org")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr", repos = "http://cran.us.r-project.org")

if (!require("msigdbr", quietly = TRUE))
  install.packages("msigdbr", repos = "http://cran.us.r-project.org")

# Bioc packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")

if (!require("missMethyl", quietly = TRUE))
  BiocManager::install("missMethyl")

if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")

if (!require("DMRcate", quietly = TRUE))
  BiocManager::install("DMRcate")

