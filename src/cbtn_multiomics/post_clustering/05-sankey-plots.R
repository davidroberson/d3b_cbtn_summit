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
# part 1 is links between molecular_subtype and multi-modal cluster
df <- anno_file %>%
  dplyr::filter(!is.na(molecular_subtype)) %>%
  dplyr::group_by(molecular_subtype, mm_cluster) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::arrange(molecular_subtype, mm_cluster, n) %>%
  dplyr::select(molecular_subtype, mm_cluster, n) %>%
  unique()

nodes <- data.frame(name = c(unique(sort(
  df$molecular_subtype
)), unique(sort(df$mm_cluster))))
nodes$df_nums <- seq(0, nrow(nodes) - 1)
links <- df %>% as.data.frame()
links$mm_cluster <- as.character(links$mm_cluster)

for (i in 1:nrow(links)) {
  links[i, "source"] <- nodes[which(nodes$name == links[i, "molecular_subtype"]), "df_nums"]
  links[i, "target"] <- nodes[which(nodes$name == links[i, "mm_cluster"]), "df_nums"]
}
links$value <- links$n
nodes$name

# part 2 is links between multi-modal cluster and CNS_region
df = anno_file %>%
  dplyr::group_by(molecular_subtype, mm_cluster, CNS_region) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::arrange(molecular_subtype, mm_cluster, CNS_region, n) %>%
  dplyr::select(molecular_subtype, mm_cluster, CNS_region, n) %>%
  unique()

nodes_2 <- data.frame(name = c(unique(sort(
  df$molecular_subtype
)), unique(sort(df$mm_cluster)), unique(sort(df$CNS_region))))
nodes_2$df_nums <- seq(0, nrow(nodes_2) - 1)
links_2 <- df %>% as.data.frame()
links_2$mm_cluster <- as.character(links_2$mm_cluster)

for (i in 1:nrow(links_2)) {
  links_2[i, "source"] <- nodes_2[which(nodes_2$name == links_2[i, "mm_cluster"]), "df_nums"]
  links_2[i, "target"] <- nodes_2[which(nodes_2$name == links_2[i, "CNS_region"]), "df_nums"]
}
links_2$value <- links_2$n
nodes_2$name

# create final set of dataframes
links_final <- rbind(
  links %>% dplyr::select(source, target, value),
  links_2 %>% dplyr::select(source, target, value)
)
nodes_final <- nodes_2

# prepare color scale: each node gets a specific color.
# remove spaces from names
nodes_final$name <- gsub(", | ", "_", nodes_final$name)
my_color <- 'd3.scaleOrdinal()
.domain(["MB_Group3", "MB_Group4", "MB_SHH", "MB_WNT",
"1", "2", "3", "4",
"5", "6", "7", "8",
"9", "10", "11", "12",
"13", "14",
"Hemispheric", "Mixed", "Posterior_fossa", "Ventricles"])
.range(["#00BFC4", "#F8766D", "#C77CFF", "#7CAE00",
"#d227da","#beaed4", "#ff1493", "#abcd21",
"#8a2be2", "#ffe135", "#228b22", "#967117",
"#21abcd", "#3a56ca", "#162e95", "#ff8c00",
"#ff4040", "#8b0000",
"#c27ba0", "#3d85c6", "#bf9000", "#b3e2cd"])'

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
saveNetwork(network = sn, file = file.path(plots_dir, "sn.html"))
webshot(
  url = file.path(plots_dir, "sn.html"),
  file = file.path(
    plots_dir,
    "molecular_subtype_vs_mm_clusters_cns_region_sankey.pdf"
  ),
  vwidth = 1000,
  vheight = 800
)

# 2) sankey plot of dkfz_v12_methylation_subclass vs multi-modal cluster
# part 1 is links between dkfz_v12_methylation_subclass and multi-modal cluster
df <- anno_file %>%
  dplyr::filter(!is.na(dkfz_v12_methylation_subclass)) %>%
  dplyr::group_by(dkfz_v12_methylation_subclass, mm_cluster) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::arrange(dkfz_v12_methylation_subclass, mm_cluster, n) %>%
  dplyr::select(dkfz_v12_methylation_subclass, mm_cluster, n) %>%
  unique()

nodes <- data.frame(name = c(unique(sort(
  df$dkfz_v12_methylation_subclass
)), unique(sort(df$mm_cluster))))
nodes$df_nums <- seq(0, nrow(nodes) - 1)
links <- df %>% as.data.frame()
links$mm_cluster <- as.character(links$mm_cluster)

for (i in 1:nrow(links)) {
  links[i, "source"] <- nodes[which(nodes$name == links[i, "dkfz_v12_methylation_subclass"]), "df_nums"]
  links[i, "target"] <- nodes[which(nodes$name == links[i, "mm_cluster"]), "df_nums"]
}
links$value <- links$n
nodes$name

# part 2 is links between multi-modal cluster and CNS_region
df <- anno_file %>%
  dplyr::group_by(dkfz_v12_methylation_subclass, mm_cluster, CNS_region) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::arrange(dkfz_v12_methylation_subclass, mm_cluster, CNS_region, n) %>%
  dplyr::select(dkfz_v12_methylation_subclass, mm_cluster, CNS_region, n) %>%
  unique()

nodes_2 <- data.frame(name = c(unique(sort(
  df$dkfz_v12_methylation_subclass
)), unique(sort(df$mm_cluster)), unique(sort(df$CNS_region))))
nodes_2$df_nums <- seq(0, nrow(nodes_2) - 1)
links_2 <- df %>% as.data.frame()
links_2$mm_cluster <- as.character(links_2$mm_cluster)

for (i in 1:nrow(links_2)) {
  links_2[i, "source"] <- nodes_2[which(nodes_2$name == links_2[i, "mm_cluster"]), "df_nums"]
  links_2[i, "target"] <- nodes_2[which(nodes_2$name == links_2[i, "CNS_region"]), "df_nums"]
}
links_2$value <- links_2$n
nodes_2$name

# create final set of dataframes
links_final <- rbind(
  links %>% dplyr::select(source, target, value),
  links_2 %>% dplyr::select(source, target, value)
)
nodes_final <- nodes_2

# prepare color scale: each node gets a specific color.
# remove spaces from names
nodes_final$name <- gsub(", | ", "_", nodes_final$name)
my_color <- 'd3.scaleOrdinal()
.domain(["ARMS", "CTRL_CBM", "ETMR_Atyp", "MB_G34_I", 
"MB_G34_II", "MB_G34_III", "MB_G34_IV", "MB_G34_V", 
"MB_G34_VI", "MB_G34_VII", "MB_G34_VIII", "MB_MYO", 
"MB_SHH_1", "MB_SHH_2", "MB_SHH_3", "MB_SHH_4", "MB_WNT",
"1", "2", "3", "4",
"5", "6", "7", "8",
"9", "10", "11", "12",
"13", "14",
"Hemispheric", "Mixed", "Posterior_fossa", "Ventricles"])
.range(["#fabed4", "#a52a2a" , "#bf9000", "#00BFC4",
"#9fc5e8", "#4363d8", "#0000ff", "#ffe135", 
"#bf9000", "#ffa500", "#F8766D", "#000000", 
"#dcbeff", "#C77CFF", "#911eb4", "#FF4DC5", "#7CAE00", 
"#d227da","#beaed4", "#ff1493", "#abcd21",
"#8a2be2", "#ffe135", "#228b22", "#967117",
"#21abcd", "#3a56ca", "#162e95", "#ff8c00",
"#ff4040", "#8b0000",
"#c27ba0", "#3d85c6", "#bf9000", "#b3e2cd"])'

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
  colourScale = my_color)

# save as an html, then as PDF
saveNetwork(network = sn, file = file.path(plots_dir, "sn.html"))
webshot(
  url = file.path(plots_dir, "sn.html"),
  file = file.path(
    plots_dir,
    "dkfz_v12_methylation_subclass_vs_mm_clusters_cns_region_sankey.pdf"
  ),
  vwidth = 1000,
  vheight = 1000
)

# remove intermediate files
system(command = paste("rm", file.path(plots_dir, "sn.html")))
system(command = paste("rm -r", file.path(plots_dir, "sn_files")))

################################### EFS event type ###################################

# 1) sankey plot of molecular_subtype vs intNMF cluster
# part 1 is links between molecular_subtype and intNMF cluster
df <- anno_file %>%
  filter(!is.na(molecular_subtype)) %>%
  group_by(molecular_subtype, mm_cluster) %>%
  mutate(n = n()) %>%
  arrange(molecular_subtype, mm_cluster, n) %>%
  dplyr::select(molecular_subtype, mm_cluster, n) %>%
  unique()

nodes <- data.frame(name = c(unique(sort(
  df$molecular_subtype
)), unique(sort(df$mm_cluster))))
nodes$df_nums <- seq(0, nrow(nodes) - 1)
links <- df %>% as.data.frame()
links$mm_cluster <- as.character(links$mm_cluster)

for (i in 1:nrow(links)) {
  links[i, "source"] <- nodes[which(nodes$name == links[i, "molecular_subtype"]), "df_nums"]
  links[i, "target"] <- nodes[which(nodes$name == links[i, "mm_cluster"]), "df_nums"]
}
links$value <- links$n
nodes$name

# part 2 is links between intNMF cluster and EFS_event_type
df = anno_file %>%
  group_by(molecular_subtype, mm_cluster, EFS_event_type) %>%
  mutate(n = n()) %>%
  arrange(molecular_subtype, mm_cluster, EFS_event_type, n) %>%
  dplyr::select(molecular_subtype, mm_cluster, EFS_event_type, n) %>%
  unique()

nodes_2 <- data.frame(name = c(unique(sort(
  df$molecular_subtype
)), unique(sort(df$mm_cluster)), unique(sort(df$EFS_event_type))))
nodes_2$df_nums <- seq(0, nrow(nodes_2) - 1)
links_2 <- df %>% as.data.frame()
links_2$mm_cluster <- as.character(links_2$mm_cluster)

for (i in 1:nrow(links_2)) {
  links_2[i, "source"] <- nodes_2[which(nodes_2$name == links_2[i, "mm_cluster"]), "df_nums"]
  links_2[i, "target"] <- nodes_2[which(nodes_2$name == links_2[i, "EFS_event_type"]), "df_nums"]
}
links_2$value <- links_2$n
nodes_2$name

# create final set of dataframes
links_final <- rbind(
  links %>% dplyr::select(source, target, value),
  links_2 %>% dplyr::select(source, target, value)
)
nodes_final <- nodes_2

# prepare color scale: each node gets a specific color.
# remove spaces from names
nodes_final$name <- gsub(", | | - ", "_", nodes_final$name)
my_color <- 'd3.scaleOrdinal()
.domain(["MB_SHH", "MB_WNT", "MB_Group3", "MB_Group4",
"1", "2", "3", "4",
"5", "6", "7", "8",
"9", "10", "11", "12",
"13", "14",
"Deceased_due_to_disease", "Not_Applicable", "Progressive", "Progressive_Metastatic",
"Recurrence", "Recurrence_Metastatic", "Second_Malignancy"]) 
.range(["#C77CFF", "#7CAE00" , "#00BFC4", "#F8766D",
"#d227da","#beaed4", "#ff1493", "#abcd21",
"#8a2be2", "#ffe135", "#228b22", "#967117",
"#21abcd", "#3a56ca", "#162e95", "#ff8c00",
"#ff4040", "#8b0000",
"#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FD7446"])'

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
saveNetwork(network = sn, file = file.path(plots_dir, "sn.html"))
webshot(
  url = file.path(plots_dir, "sn.html"),
  file = file.path(
    plots_dir,
    "molecular_subtype_vs_mm_clusters_efs_event_sankey.pdf"
  ),
  vwidth = 1000,
  vheight = 800
)

# 2) sankey plot of dkfz_v12_methylation_subclass vs intNMF cluster
# part 1 is links between dkfz_v12_methylation_subclass and intNMF cluster
df <- anno_file %>%
  dplyr::filter(!is.na(dkfz_v12_methylation_subclass)) %>%
  dplyr::group_by(dkfz_v12_methylation_subclass, mm_cluster) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::arrange(dkfz_v12_methylation_subclass, mm_cluster, n) %>%
  dplyr::select(dkfz_v12_methylation_subclass, mm_cluster, n) %>%
  unique()

nodes <- data.frame(name = c(unique(sort(
  df$dkfz_v12_methylation_subclass
)), unique(sort(df$mm_cluster))))
nodes$df_nums <- seq(0, nrow(nodes) - 1)
links <- df %>% as.data.frame()
links$mm_cluster <- as.character(links$mm_cluster)

for (i in 1:nrow(links)) {
  links[i, "source"] <- nodes[which(nodes$name == links[i, "dkfz_v12_methylation_subclass"]), "df_nums"]
  links[i, "target"] <- nodes[which(nodes$name == links[i, "mm_cluster"]), "df_nums"]
}
links$value <- links$n
nodes$name

# part 2 is links between intNMF cluster and EFS_event_type
df <- anno_file %>%
  dplyr::group_by(dkfz_v12_methylation_subclass, mm_cluster, EFS_event_type) %>%
  dplyr::mutate(n = n()) %>%
  dplyr::arrange(dkfz_v12_methylation_subclass, mm_cluster, EFS_event_type, n) %>%
  dplyr::select(dkfz_v12_methylation_subclass, mm_cluster, EFS_event_type, n) %>%
  unique()

nodes_2 <- data.frame(name = c(unique(sort(
  df$dkfz_v12_methylation_subclass
)), unique(sort(df$mm_cluster)), unique(sort(df$EFS_event_type))))
nodes_2$df_nums <- seq(0, nrow(nodes_2) - 1)
links_2 <- df %>% as.data.frame()
links_2$mm_cluster <- as.character(links_2$mm_cluster)

for (i in 1:nrow(links_2)) {
  links_2[i, "source"] <- nodes_2[which(nodes_2$name == links_2[i, "mm_cluster"]), "df_nums"]
  links_2[i, "target"] <- nodes_2[which(nodes_2$name == links_2[i, "EFS_event_type"]), "df_nums"]
}
links_2$value <- links_2$n
nodes_2$name

# create final set of dataframes
links_final <- rbind(
  links %>% dplyr::select(source, target, value),
  links_2 %>% dplyr::select(source, target, value)
)
nodes_final <- nodes_2

# prepare color scale: each node gets a specific color.
# remove spaces from names
nodes_final$name <- gsub(", | | - ", "_", nodes_final$name)
my_color <- 'd3.scaleOrdinal()
.domain(["ARMS", "CTRL_CBM", "ETMR_Atyp", "MB_G34_I", 
"MB_G34_II", "MB_G34_III", "MB_G34_IV", "MB_G34_V", 
"MB_G34_VI", "MB_G34_VII", "MB_G34_VIII", "MB_MYO", 
"MB_SHH_1", "MB_SHH_2", "MB_SHH_3", "MB_SHH_4", "MB_WNT",
"1", "2", "3", "4",
"5", "6", "7", "8",
"9", "10", "11", "12",
"13", "14",
"Deceased_due_to_disease", "Not_Applicable", "Progressive", "Progressive_Metastatic",
"Recurrence", "Recurrence_Metastatic", "Second_Malignancy"])
.range(["#fabed4", "#a52a2a" , "#bf9000", "#00BFC4",
"#9fc5e8", "#4363d8", "#0000ff", "#ffe135", 
"#bf9000", "#ffa500", "#F8766D", "#000000", 
"#dcbeff", "#C77CFF", "#911eb4", "#FF4DC5", "#7CAE00", 
"#d227da","#beaed4", "#ff1493", "#abcd21",
"#8a2be2", "#ffe135", "#228b22", "#967117",
"#21abcd", "#3a56ca", "#162e95", "#ff8c00",
"#ff4040", "#8b0000",
"#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FD7446"])'

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
saveNetwork(network = sn, file = file.path(plots_dir, "sn.html"))
webshot(
  url = file.path(plots_dir, "sn.html"),
  file = file.path(
    plots_dir,
    "dkfz_v12_methylation_subclass_vs_mm_clusters_efs_event_sankey.pdf"
  ),
  vwidth = 1000,
  vheight = 1000
)

# remove intermediate files
system(command = paste("rm", file.path(plots_dir, "sn.html")))
system(command = paste("rm -r", file.path(plots_dir, "sn_files")))
