# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", repos = "http://cran.us.r-project.org")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")

if (!require("datawizard", quietly = TRUE))
  install.packages("datawizard", repos = "http://cran.us.r-project.org")

if (!require("reshape2", quietly = TRUE))
  install.packages("reshape2", repos = "http://cran.us.r-project.org")

# Bioc packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")

if (!require("rtracklayer", quietly = TRUE))
  BiocManager::install("rtracklayer")
