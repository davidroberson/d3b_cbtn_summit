# script to install required packages for the module

# CRAN packages
if (!require("optparse", quietly = TRUE))
  install.packages("optparse")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!require("corrplot", quietly = TRUE))
  install.packages("corrplot")

if (!require("data.table", quietly = TRUE))
  install.packages("data.table")

if (!require("ggplotify", quietly = TRUE))
  install.packages("ggplotify")

if (!require("ggpubr", quietly = TRUE))
  install.packages("ggpubr")

if (!require("gplots", quietly = TRUE))
  install.packages("gplots")

if (!require("gridExtra", quietly = TRUE))
  install.packages("gridExtra")

if (!require("gtsummary", quietly = TRUE))
  install.packages("gtsummary")

if (!require("mclust", quietly = TRUE))
  install.packages("mclust")

if (!require("networkD3", quietly = TRUE))
  install.packages("networkD3")

if (!require("pheatmap", quietly = TRUE))
  install.packages("pheatmap")

if (!require("randomcoloR", quietly = TRUE))
  install.packages("randomcoloR")

if (!require("survival", quietly = TRUE))
  install.packages("survival")

if (!require("survminer", quietly = TRUE))
  install.packages("survminer")

if (!require("webshot", quietly = TRUE))
  install.packages("webshot")
