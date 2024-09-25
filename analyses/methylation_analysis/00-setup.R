# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!require("data.table", quietly = TRUE))
  install.packages("data.table")

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!require("msigdbr", quietly = TRUE))
  install.packages("msigdbr")

# Bioc packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("missMethyl", quietly = TRUE))
  BiocManager::install("missMethyl")

if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")

if (!require("DMRcate", quietly = TRUE))
  BiocManager::install("DMRcate")

