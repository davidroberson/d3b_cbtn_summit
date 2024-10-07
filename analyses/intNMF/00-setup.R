# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", repos = "http://cran.us.r-project.org")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")

if (!require("IntNMF", quietly = TRUE))
  install.packages("IntNMF", repos = "http://cran.us.r-project.org")
