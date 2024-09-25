# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!require("IntNMF", quietly = TRUE))
  install.packages("IntNMF")
