# script to generate sankey plots
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(networkD3)
  library(webshot)
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

# combine multi-modal clusters with RNA-derived molecular subtypes
histology_file <- opt$histology_file
anno_file_rna <- read_tsv(file = histology_file) %>%
  dplyr::select(Kids_First_Biospecimen_ID,
                molecular_subtype,
                CNS_region,
                EFS_event_type) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID") %>%
  unique() %>%
  dplyr::inner_join(mm_clusters, by = "Kids_First_Biospecimen_ID_RNA")

# combine multi-modal clusters with methylation-derived subclass
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

############################# generate sankey plots ############################
################################### CNS_region ###################################

# 1) sankey plot of molecular_subtype vs multi-modal cluster
# Create a generic function to create a placeholder PDF with a message
create_placeholder_pdf <- function(file_path, message) {
  pdf(file = file_path, width = 8, height = 6)
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = message)
  text(1, 1, message, cex = 1.5)
  dev.off()
}

# Function to safely create a sankey plot
create_sankey_plot <- function(part1_data, part2_data, part1_src_col, part1_tgt_col, part2_src_col, part2_tgt_col, 
                               output_file, title, domain_vals = NULL, range_vals = NULL) {
  tryCatch({
    # Part 1 processing
    df1 <- part1_data %>%
      dplyr::filter(!is.na(!!rlang::sym(part1_src_col))) %>%
      dplyr::group_by(!!rlang::sym(part1_src_col), !!rlang::sym(part1_tgt_col)) %>%
      dplyr::mutate(n = n()) %>%
      dplyr::arrange(!!rlang::sym(part1_src_col), !!rlang::sym(part1_tgt_col), n) %>%
      dplyr::select(!!rlang::sym(part1_src_col), !!rlang::sym(part1_tgt_col), n) %>%
      unique()
    
    # Check if we have data
    if (nrow(df1) == 0) {
      cat("Warning: No data available for part 1 of sankey plot", title, "\n")
      create_placeholder_pdf(output_file, paste("No data available for part 1 of sankey plot", title))
      return(FALSE)
    }
    
    # Create nodes and links for part 1
    src_vals <- unique(sort(df1[[part1_src_col]]))
    tgt_vals <- unique(sort(df1[[part1_tgt_col]]))
    
    # Convert to character
    df1[[part1_tgt_col]] <- as.character(df1[[part1_tgt_col]])
    
    # Create nodes dataframe
    nodes <- data.frame(name = c(src_vals, tgt_vals))
    nodes$df_nums <- seq(0, nrow(nodes) - 1)
    
    # Create links dataframe
    links <- df1 %>% as.data.frame()
    links$source <- NA
    links$target <- NA
    
    # Safely assign source and target indices
    for (i in 1:nrow(links)) {
      src_idx <- which(nodes$name == links[i, part1_src_col])
      tgt_idx <- which(nodes$name == links[i, part1_tgt_col])
      
      # Check if indices were found
      if (length(src_idx) > 0 && length(tgt_idx) > 0) {
        links[i, "source"] <- nodes$df_nums[src_idx[1]]
        links[i, "target"] <- nodes$df_nums[tgt_idx[1]]
      } else {
        cat("Warning: Could not find node for row", i, "in part 1\n")
        if (length(src_idx) == 0) cat("  Missing source node:", links[i, part1_src_col], "\n")
        if (length(tgt_idx) == 0) cat("  Missing target node:", links[i, part1_tgt_col], "\n")
      }
    }
    links$value <- links$n
    
    # Part 2 processing
    df2 <- part2_data %>%
      dplyr::group_by(!!rlang::sym(part1_src_col), !!rlang::sym(part1_tgt_col), !!rlang::sym(part2_tgt_col)) %>%
      dplyr::mutate(n = n()) %>%
      dplyr::arrange(!!rlang::sym(part1_src_col), !!rlang::sym(part1_tgt_col), !!rlang::sym(part2_tgt_col), n) %>%
      dplyr::select(!!rlang::sym(part1_src_col), !!rlang::sym(part1_tgt_col), !!rlang::sym(part2_tgt_col), n) %>%
      unique()
    
    # Check if we have data
    if (nrow(df2) == 0) {
      cat("Warning: No data available for part 2 of sankey plot", title, "\n")
      create_placeholder_pdf(output_file, paste("No data available for part 2 of sankey plot", title))
      return(FALSE)
    }
    
    # Create nodes and links for part 2
    src_vals <- unique(sort(df2[[part1_src_col]]))
    mid_vals <- unique(sort(df2[[part1_tgt_col]]))
    tgt_vals <- unique(sort(df2[[part2_tgt_col]]))
    
    # Convert to character
    df2[[part1_tgt_col]] <- as.character(df2[[part1_tgt_col]])
    
    # Create nodes dataframe
    nodes2 <- data.frame(name = c(src_vals, mid_vals, tgt_vals))
    nodes2$df_nums <- seq(0, nrow(nodes2) - 1)
    
    # Create links dataframe
    links2 <- df2 %>% as.data.frame()
    links2$source <- NA
    links2$target <- NA
    
    # Safely assign source and target indices
    for (i in 1:nrow(links2)) {
      src_idx <- which(nodes2$name == links2[i, part2_src_col])
      tgt_idx <- which(nodes2$name == links2[i, part2_tgt_col])
      
      # Check if indices were found
      if (length(src_idx) > 0 && length(tgt_idx) > 0) {
        links2[i, "source"] <- nodes2$df_nums[src_idx[1]]
        links2[i, "target"] <- nodes2$df_nums[tgt_idx[1]]
      } else {
        cat("Warning: Could not find node for row", i, "in part 2\n")
        if (length(src_idx) == 0) cat("  Missing source node:", links2[i, part2_src_col], "\n")
        if (length(tgt_idx) == 0) cat("  Missing target node:", links2[i, part2_tgt_col], "\n")
      }
    }
    links2$value <- links2$n
    
    # Remove NA links
    links <- links[!is.na(links$source) & !is.na(links$target),]
    links2 <- links2[!is.na(links2$source) & !is.na(links2$target),]
    
    # Check if we still have data
    if (nrow(links) == 0 || nrow(links2) == 0) {
      cat("Warning: No valid links available for sankey plot", title, "\n")
      create_placeholder_pdf(output_file, paste("No valid links available for sankey plot", title))
      return(FALSE)
    }
    
    # create final set of dataframes
    links_final <- rbind(
      links %>% dplyr::select(source, target, value),
      links2 %>% dplyr::select(source, target, value)
    )
    nodes_final <- nodes2
    
    # prepare color scale: each node gets a specific color.
    # remove spaces from names
    nodes_final$name <- gsub(", | ", "_", nodes_final$name)
    
    # Use provided domain values or defaults
    if (is.null(domain_vals)) {
      domain_vals <- paste0('"', paste(nodes_final$name, collapse = '", "'), '"')
    } else {
      # Ensure all node names are in the domain
      all_names <- unique(nodes_final$name)
      missing_names <- all_names[!all_names %in% domain_vals]
      if (length(missing_names) > 0) {
        cat("Warning: Some node names are not in the domain:", paste(missing_names, collapse = ", "), "\n")
        domain_vals <- c(domain_vals, missing_names)
      }
      domain_vals <- paste0('"', paste(domain_vals, collapse = '", "'), '"')
    }
    
    # Use provided range values or defaults
    if (is.null(range_vals)) {
      # Generate a reasonable number of colors
      n_colors <- length(unique(nodes_final$name))
      range_vals <- paste0('"', paste(rainbow(n_colors), collapse = '", "'), '"')
    } else {
      range_vals <- paste0('"', paste(range_vals, collapse = '", "'), '"')
    }
    
    my_color <- paste0('d3.scaleOrdinal()\n.domain([', domain_vals, '])\n.range([', range_vals, '])')
    
    # plot
    sn <- sankeyNetwork(
      Links = links_final,
      Nodes = nodes_final,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      units = "TWh",
      fontSize = 16,
      nodeWidth = 30,
      colourScale = my_color
    )
    
    # save as an html, then as PDF
    temp_html <- file.path(plots_dir, "sn_temp.html")
    
    tryCatch({
      saveNetwork(network = sn, file = temp_html)
      webshot(
        url = temp_html,
        file = output_file,
        vwidth = 1000,
        vheight = 800
      )
      # Clean up temporary files
      if (file.exists(temp_html)) file.remove(temp_html)
      temp_dir <- file.path(plots_dir, "sn_temp_files")
      if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
      return(TRUE)
    }, error = function(e) {
      cat("Error saving sankey plot:", e$message, "\n")
      create_placeholder_pdf(output_file, paste("Error creating sankey plot:", e$message))
      return(FALSE)
    })
    
  }, error = function(e) {
    cat("Error creating sankey plot:", e$message, "\n")
    create_placeholder_pdf(output_file, paste("Error creating sankey plot:", e$message))
    return(FALSE)
  })
}

# Prepare domain and range values for color scales
molecular_domain <- c("MB_Group3", "MB_Group4", "MB_SHH", "MB_WNT",
                      "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
                      "11", "12", "13", "14",
                      "Hemispheric", "Mixed", "Posterior_fossa", "Ventricles")

molecular_range <- c("#00BFC4", "#F8766D", "#C77CFF", "#7CAE00",
                     "#d227da","#beaed4", "#ff1493", "#abcd21",
                     "#8a2be2", "#ffe135", "#228b22", "#967117",
                     "#21abcd", "#3a56ca", "#162e95", "#ff8c00",
                     "#ff4040", "#8b0000",
                     "#c27ba0", "#3d85c6", "#bf9000", "#b3e2cd")

# Create the first sankey plot
cat("Creating molecular subtype vs mm_clusters sankey plot\n")
create_sankey_plot(
  part1_data = anno_file, 
  part2_data = anno_file, 
  part1_src_col = "molecular_subtype", 
  part1_tgt_col = "mm_cluster", 
  part2_src_col = "mm_cluster", 
  part2_tgt_col = "CNS_region", 
  output_file = file.path(plots_dir, "molecular_subtype_vs_mm_clusters_cns_region_sankey.pdf"),
  title = "Molecular subtype vs mm_clusters",
  domain_vals = molecular_domain,
  range_vals = molecular_range
)

# 2) sankey plot of dkfz_v12_methylation_subclass vs multi-modal cluster
# Prepare domain and range values for dkfz color scales
dkfz_domain <- c("ARMS", "CTRL_CBM", "ETMR_Atyp", "MB_G34_I", 
              "MB_G34_II", "MB_G34_III", "MB_G34_IV", "MB_G34_V", 
              "MB_G34_VI", "MB_G34_VII", "MB_G34_VIII", "MB_MYO", 
              "MB_SHH_1", "MB_SHH_2", "MB_SHH_3", "MB_SHH_4", "MB_WNT",
              "1", "2", "3", "4", "5", "6", "7", "8",
              "9", "10", "11", "12", "13", "14",
              "Hemispheric", "Mixed", "Posterior_fossa", "Ventricles")

dkfz_range <- c("#fabed4", "#a52a2a" , "#bf9000", "#00BFC4",
              "#9fc5e8", "#4363d8", "#0000ff", "#ffe135", 
              "#bf9000", "#ffa500", "#F8766D", "#000000", 
              "#dcbeff", "#C77CFF", "#911eb4", "#FF4DC5", "#7CAE00", 
              "#d227da","#beaed4", "#ff1493", "#abcd21",
              "#8a2be2", "#ffe135", "#228b22", "#967117",
              "#21abcd", "#3a56ca", "#162e95", "#ff8c00",
              "#ff4040", "#8b0000",
              "#c27ba0", "#3d85c6", "#bf9000", "#b3e2cd")

# Create the second sankey plot
cat("Creating dkfz_v12_methylation_subclass vs mm_clusters sankey plot\n")
create_sankey_plot(
  part1_data = anno_file, 
  part2_data = anno_file, 
  part1_src_col = "dkfz_v12_methylation_subclass", 
  part1_tgt_col = "mm_cluster", 
  part2_src_col = "mm_cluster", 
  part2_tgt_col = "CNS_region", 
  output_file = file.path(plots_dir, "dkfz_v12_methylation_subclass_vs_mm_clusters_cns_region_sankey.pdf"),
  title = "DKFZ v12 methylation subclass vs mm_clusters",
  domain_vals = dkfz_domain,
  range_vals = dkfz_range
)

# We'll safely remove temporary files created by our function
# Our function already handles this internally
# but let's ensure any leftover original files are cleaned up
tryCatch({
  # Try to remove any old files that might exist
  if (file.exists(file.path(plots_dir, "sn.html"))) {
    file.remove(file.path(plots_dir, "sn.html"))
  }
  
  sn_files_dir <- file.path(plots_dir, "sn_files")
  if (dir.exists(sn_files_dir)) {
    unlink(sn_files_dir, recursive = TRUE)
  }
}, error = function(e) {
  cat("Warning: Failed to clean up temporary files:", e$message, "\n")
})

################################### EFS event type ###################################

# Prepare domain and range values for EFS event type sankey plots
efs_molecular_domain <- c("MB_SHH", "MB_WNT", "MB_Group3", "MB_Group4",
                          "1", "2", "3", "4", "5", "6", "7", "8",
                          "9", "10", "11", "12", "13", "14",
                          "Deceased_due_to_disease", "Not_Applicable", "Progressive", 
                          "Progressive_Metastatic", "Recurrence", "Recurrence_Metastatic", 
                          "Second_Malignancy")

efs_molecular_range <- c("#C77CFF", "#7CAE00" , "#00BFC4", "#F8766D",
                       "#d227da","#beaed4", "#ff1493", "#abcd21",
                       "#8a2be2", "#ffe135", "#228b22", "#967117",
                       "#21abcd", "#3a56ca", "#162e95", "#ff8c00",
                       "#ff4040", "#8b0000",
                       "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FD7446")

# Create the third sankey plot
cat("Creating molecular subtype vs mm_clusters for EFS event type sankey plot\n")
create_sankey_plot(
  part1_data = anno_file, 
  part2_data = anno_file, 
  part1_src_col = "molecular_subtype", 
  part1_tgt_col = "mm_cluster", 
  part2_src_col = "mm_cluster", 
  part2_tgt_col = "EFS_event_type", 
  output_file = file.path(plots_dir, "molecular_subtype_vs_mm_clusters_efs_event_sankey.pdf"),
  title = "Molecular subtype vs mm_clusters (EFS event)",
  domain_vals = efs_molecular_domain,
  range_vals = efs_molecular_range
)

# Prepare domain and range values for DKFZ + EFS event type
efs_dkfz_domain <- c("ARMS", "CTRL_CBM", "ETMR_Atyp", "MB_G34_I", 
                   "MB_G34_II", "MB_G34_III", "MB_G34_IV", "MB_G34_V", 
                   "MB_G34_VI", "MB_G34_VII", "MB_G34_VIII", "MB_MYO", 
                   "MB_SHH_1", "MB_SHH_2", "MB_SHH_3", "MB_SHH_4", "MB_WNT",
                   "1", "2", "3", "4", "5", "6", "7", "8",
                   "9", "10", "11", "12", "13", "14",
                   "Deceased_due_to_disease", "Not_Applicable", "Progressive", 
                   "Progressive_Metastatic", "Recurrence", "Recurrence_Metastatic", 
                   "Second_Malignancy")

efs_dkfz_range <- c("#fabed4", "#a52a2a" , "#bf9000", "#00BFC4",
                  "#9fc5e8", "#4363d8", "#0000ff", "#ffe135", 
                  "#bf9000", "#ffa500", "#F8766D", "#000000", 
                  "#dcbeff", "#C77CFF", "#911eb4", "#FF4DC5", "#7CAE00", 
                  "#d227da","#beaed4", "#ff1493", "#abcd21",
                  "#8a2be2", "#ffe135", "#228b22", "#967117",
                  "#21abcd", "#3a56ca", "#162e95", "#ff8c00",
                  "#ff4040", "#8b0000",
                  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FD7446")

# Create the fourth sankey plot
cat("Creating dkfz_v12_methylation_subclass vs mm_clusters for EFS event type sankey plot\n")
create_sankey_plot(
  part1_data = anno_file, 
  part2_data = anno_file, 
  part1_src_col = "dkfz_v12_methylation_subclass", 
  part1_tgt_col = "mm_cluster", 
  part2_src_col = "mm_cluster", 
  part2_tgt_col = "EFS_event_type", 
  output_file = file.path(plots_dir, "dkfz_v12_methylation_subclass_vs_mm_clusters_efs_event_sankey.pdf"),
  title = "DKFZ v12 methylation subclass vs mm_clusters (EFS event)",
  domain_vals = efs_dkfz_domain,
  range_vals = efs_dkfz_range
)

# Final cleanup of any temporary files
tryCatch({
  # Try to remove any files that might exist
  if (file.exists(file.path(plots_dir, "sn.html"))) {
    file.remove(file.path(plots_dir, "sn.html"))
  }
  
  if (file.exists(file.path(plots_dir, "sn_temp.html"))) {
    file.remove(file.path(plots_dir, "sn_temp.html"))
  }
  
  sn_files_dir <- file.path(plots_dir, "sn_files")
  if (dir.exists(sn_files_dir)) {
    unlink(sn_files_dir, recursive = TRUE)
  }
  
  sn_temp_files_dir <- file.path(plots_dir, "sn_temp_files")
  if (dir.exists(sn_temp_files_dir)) {
    unlink(sn_temp_files_dir, recursive = TRUE)
  }
}, error = function(e) {
  cat("Warning: Failed to clean up temporary files:", e$message, "\n")
})
