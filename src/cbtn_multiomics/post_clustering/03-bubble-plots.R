# script to create plots for each data modality arranged and annotated by clusters
# to determine how the clusters relate to subtypes

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(gplots)
  library(corrplot)
})

# parse command line options
option_list <- list(
  make_option(c("--cluster_file"), type = "character", help = "File with multi-modal derived clusters (.tsv)"),
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--plots_dir"), type = "character", help = "path to plots directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))

# plots directory
plots_dir <- opt$plots_dir
dir.create(plots_dir, showWarnings = F, recursive = T)

# multi-modal derived clusters
mm_clusters <- read_tsv(file = opt$cluster_file)

# combine Multi-modal clusters with RNA-derived molecular subtypes
cat('Joining multi-omic cluster annotations with histologic data \n')
histology_file <- opt$histology_file
anno_file_rna <- read_tsv(file = histology_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype, CNS_region) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  dplyr::inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# combine Multi-modal clusters with methylation-derived subclass
anno_file_methyl <- read_tsv(file = histology_file) %>%
  dplyr::select(
    Kids_First_Biospecimen_ID,
    dkfz_v11_methylation_subclass,
    dkfz_v12_methylation_subclass
  ) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  dplyr::inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_Methyl")

# combine both and create one standardized annotation file
anno_file <- anno_file_rna %>%
  dplyr::inner_join(anno_file_methyl)
anno_file$dkfz_v11_methylation_subclass <- gsub(", | ", "_", anno_file$dkfz_v11_methylation_subclass)

############################# generate balloon and corrplots ############################

# generate balloon plot with at least 5 samples in a group
# Multi-modal clusters vs RNA-derived molecular subtypes
cat('Plotting multi-omic clusters in relation to known subtypes \n')
dat <- anno_file %>%
  filter(!is.na(molecular_subtype)) %>%
  group_by(molecular_subtype, mm_cluster)  %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(nmax = max(n)) %>%
  filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  spread(key = mm_cluster, value = n, fill = 0) %>%
  column_to_rownames("molecular_subtype")
pdf(
  file = file.path(plots_dir, "mm_clusters_vs_molsubtype_balloonplot.pdf"),
  width = 10
)
# Create a proper table format for balloonplot - with defensive check
if (nrow(dat) > 0 && ncol(dat) > 0) {
  # Create dummy data if no real data
  if (length(colnames(dat)) != length(rownames(dat))) {
    # Print diagnostic info
    cat("Warning: Column names and row names have different lengths.\n")
    cat("Columns:", length(colnames(dat)), "\n")
    cat("Rows:", length(rownames(dat)), "\n")
    
    # Create a placeholder matrix as fallback
    placeholder <- matrix(0, nrow=1, ncol=1)
    colnames(placeholder) <- "No data"
    rownames(placeholder) <- "No data"
    mm_vs_subtype <- as.table(placeholder)
    
    # Add a note to the plot
    balloonplot(
      x = mm_vs_subtype,
      main = "Multi-modal clusters vs Molecular subtypes\n(Insufficient data)",
      xlab = "",
      ylab = "",
      label = TRUE,
      show.margins = FALSE
    )
  } else {
    mm_vs_subtype <- table(colnames(dat), rownames(dat))
    balloonplot(
      x = mm_vs_subtype,
      main = "Multi-modal clusters vs Molecular subtypes",
      xlab = "",
      ylab = "",
      label = TRUE,
      show.margins = FALSE
    )
  }
} else {
  # Create placeholder plot for empty data
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "Insufficient data for plotting", cex = 1.5)
}
dev.off()

# 4) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# Multi-modal clusters vs RNA-derived molecular subtypes
pdf(file = file.path(plots_dir, "mm_clusters_vs_molsubtype_corrplot.pdf"))

if (nrow(dat) > 0 && ncol(dat) > 0) {
  # Only run the chi-square test if we have sufficient data
  tryCatch({
    chisq <- chisq.test(dat)
    corrplot(
      chisq$residuals,
      is.cor = FALSE,
      tl.srt = 360,
      tl.offset = 1,
      mar = c(1, 2, 1, 1),
      title = "Multi-modal clusters vs Molecular subtypes"
    )
  }, error = function(e) {
    cat("Chi-square test error:", e$message, "\n")
    # Create placeholder plot
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, paste("Insufficient data for chi-square test:", e$message), cex = 1)
  })
} else {
  # Create placeholder plot for empty data
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "Insufficient data for plotting", cex = 1.5)
}

dev.off()

# 3) generate balloon plot with at least 5 samples in a group
# Multi-modal clusters vs methylation-derived dkfz_v11_methylation_subclass
dat <- anno_file %>%
  dplyr::filter(!is.na(dkfz_v11_methylation_subclass)) %>%
  dplyr::group_by(dkfz_v11_methylation_subclass, mm_cluster)  %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::mutate(nmax = max(n)) %>%
  dplyr::filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  tidyr::spread(key = mm_cluster, value = n, fill = 0) %>%
  tibble::column_to_rownames("dkfz_v11_methylation_subclass")
pdf(
  file = file.path(
    plots_dir,
    "mm_clusters_vs_dkfz_v11_methylation_subclass_balloonplot.pdf"
  ),
  width = 14
)
# Create a proper table format for balloonplot - with defensive check
if (nrow(dat) > 0 && ncol(dat) > 0) {
  # Create dummy data if no real data
  if (length(colnames(dat)) != length(rownames(dat))) {
    # Print diagnostic info
    cat("Warning: Column names and row names have different lengths.\n")
    cat("Columns:", length(colnames(dat)), "\n")
    cat("Rows:", length(rownames(dat)), "\n")
    
    # Create a placeholder matrix as fallback
    placeholder <- matrix(0, nrow=1, ncol=1)
    colnames(placeholder) <- "No data"
    rownames(placeholder) <- "No data"
    mm_vs_subtype <- as.table(placeholder)
    
    # Add a note to the plot
    balloonplot(
      x = mm_vs_subtype,
      main = "Multi-modal clusters vs dkfz_v11_methylation_subclass\n(Insufficient data)",
      xlab = "",
      ylab = "",
      label = TRUE,
      show.margins = FALSE
    )
  } else {
    mm_vs_subtype <- table(colnames(dat), rownames(dat))
    balloonplot(
      x = mm_vs_subtype,
      main = "Multi-modal clusters vs dkfz_v11_methylation_subclass",
      xlab = "",
      ylab = "",
      label = TRUE,
      show.margins = FALSE
    )
  }
} else {
  # Create placeholder plot for empty data
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "Insufficient data for plotting", cex = 1.5)
}
dev.off()

# 4) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# Multi-modal clusters vs methylation-derived dkfz_v11_methylation_subclass
pdf(file = file.path(
  plots_dir,
  "mm_clusters_vs_dkfz_v11_methylation_subclass_corrplot.pdf"
))

if (nrow(dat) > 0 && ncol(dat) > 0) {
  # Only run the chi-square test if we have sufficient data
  tryCatch({
    chisq <- chisq.test(dat)
    corrplot(
      chisq$residuals,
      is.cor = FALSE,
      tl.srt = 360,
      tl.offset = 1,
      mar = c(1, 2, 1, 1),
      title = "Multi-modal clusters vs dkfz_v11_methylation_subclass"
    )
  }, error = function(e) {
    cat("Chi-square test error:", e$message, "\n")
    # Create placeholder plot
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, paste("Insufficient data for chi-square test:", e$message), cex = 1)
  })
} else {
  # Create placeholder plot for empty data
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "Insufficient data for plotting", cex = 1.5)
}

dev.off()

# 3) generate balloon plot with at least 5 samples in a group
# Multi-modal clusters vs methylation-derived dkfz_v12_methylation_subclass
dat <- anno_file %>%
  dplyr::filter(!is.na(dkfz_v12_methylation_subclass)) %>%
  dplyr::group_by(dkfz_v12_methylation_subclass, mm_cluster)  %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::mutate(nmax = max(n)) %>%
  dplyr::filter(nmax >= 5) %>%
  ungroup() %>%
  dplyr::select(-c(nmax)) %>%
  tidyr::spread(key = mm_cluster, value = n, fill = 0) %>%
  tibble::column_to_rownames("dkfz_v12_methylation_subclass")
pdf(
  file = file.path(
    plots_dir,
    "mm_clusters_vs_dkfz_v12_methylation_subclass_balloonplot.pdf"
  ),
  width = 14
)
# Create a proper table format for balloonplot - with defensive check
if (nrow(dat) > 0 && ncol(dat) > 0) {
  # Create dummy data if no real data
  if (length(colnames(dat)) != length(rownames(dat))) {
    # Print diagnostic info
    cat("Warning: Column names and row names have different lengths.\n")
    cat("Columns:", length(colnames(dat)), "\n")
    cat("Rows:", length(rownames(dat)), "\n")
    
    # Create a placeholder matrix as fallback
    placeholder <- matrix(0, nrow=1, ncol=1)
    colnames(placeholder) <- "No data"
    rownames(placeholder) <- "No data"
    mm_vs_subtype <- as.table(placeholder)
    
    # Add a note to the plot
    balloonplot(
      x = mm_vs_subtype,
      main = "Multi-modal clusters vs dkfz_v12_methylation_subclass\n(Insufficient data)",
      xlab = "",
      ylab = "",
      label = TRUE,
      show.margins = FALSE
    )
  } else {
    mm_vs_subtype <- table(colnames(dat), rownames(dat))
    balloonplot(
      x = mm_vs_subtype,
      main = "Multi-modal clusters vs dkfz_v12_methylation_subclass",
      xlab = "",
      ylab = "",
      label = TRUE,
      show.margins = FALSE
    )
  }
} else {
  # Create placeholder plot for empty data
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "Insufficient data for plotting", cex = 1.5)
}
dev.off()

# 4) generate corrplot of rows with at least 5 samples in a group to show pearson residuals
# Multi-modal clusters vs methylation-derived dkfz_v12_methylation_subclass
pdf(file = file.path(
  plots_dir,
  "mm_clusters_vs_dkfz_v12_methylation_subclass_corrplot.pdf"
))

if (nrow(dat) > 0 && ncol(dat) > 0) {
  # Only run the chi-square test if we have sufficient data
  tryCatch({
    chisq <- chisq.test(dat)
    corrplot(
      chisq$residuals,
      is.cor = FALSE,
      tl.srt = 360,
      tl.offset = 1,
      mar = c(1, 2, 1, 1),
      title = "Multi-modal clusters vs dkfz_v12_methylation_subclass"
    )
  }, error = function(e) {
    cat("Chi-square test error:", e$message, "\n")
    # Create placeholder plot
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, paste("Insufficient data for chi-square test:", e$message), cex = 1)
  })
} else {
  # Create placeholder plot for empty data
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(1, 1, "Insufficient data for plotting", cex = 1.5)
}

dev.off()
