# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse", repos = "http://cran.us.r-project.org")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")

if (!require("corrplot", quietly = TRUE))
  install.packages("corrplot", repos = "http://cran.us.r-project.org")

if (!require("data.table", quietly = TRUE))
  install.packages("data.table", repos = "http://cran.us.r-project.org")

if (!require("ggplotify", quietly = TRUE))
  install.packages("ggplotify", repos = "http://cran.us.r-project.org")

if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr", repos = "http://cran.us.r-project.org")

if (!require("gplots", quietly = TRUE))
  install.packages("gplots", repos = "http://cran.us.r-project.org")

if (!require("gridExtra", quietly = TRUE))
  install.packages("gridExtra", repos = "http://cran.us.r-project.org")

if (!require("gtsummary", quietly = TRUE))
  install.packages("gtsummary", repos = "http://cran.us.r-project.org")

if (!require("mclust", quietly = TRUE))
  install.packages("mclust", repos = "http://cran.us.r-project.org")

if (!require("networkD3", quietly = TRUE))
  install.packages("networkD3", repos = "http://cran.us.r-project.org")

if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap", repos = "http://cran.us.r-project.org")

if (!require("randomcoloR", quietly = TRUE))
  install.packages("randomcoloR", repos = "http://cran.us.r-project.org")

if (!require("survival", quietly = TRUE))
  install.packages("survival", repos = "http://cran.us.r-project.org")

if (!require("survminer", quietly = TRUE))
  install.packages("survminer", repos = "http://cran.us.r-project.org")

if (!require("webshot", quietly = TRUE))
  install.packages("webshot", repos = "http://cran.us.r-project.org")
