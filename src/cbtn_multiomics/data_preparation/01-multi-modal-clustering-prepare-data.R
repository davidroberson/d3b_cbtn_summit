# prepare files for multi-modal clustering
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(datawizard)
  library(reshape2)
  library(rtracklayer)
})

mem.maxVSize(vsize = 102400)

# parse command line options
option_list <- list(
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--short_histology"), type = "character", help = "Short histology of interest"),
  make_option(c("--count_file"), type = "character", help = "Expression matrix, preferably counts (.rds)"),
  make_option(c("--methyl_file"), type = "character", help = "Methylation values, preferably beta values (.rds)"),
  make_option(c("--splice_file"), type = "character", help = "PSI values (.rds)"),
  make_option(c("--gtf_file"), type = "character", help = "Gencode gtf file"),
  make_option(c("--num_features"), type = "integer", default = 1000, help = "Number of top variable features to select"),
  make_option(c("--output_dir"), type = "character", help = "path to results directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
short_histology_of_interest <- opt$short_histology

# output directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# Define the filter_cnv function inline since we can't rely on working directory
filter_cnv <- function(cnv_mat) {
  # Filter out high CNV segments
  # This is a simplified version for testing
  return(cnv_mat)
}

# read histology file and filter to short histology of interest
cat('Reading histology file \n')
histology_file <- opt$histology_file
histology_file <- readr::read_tsv(file = histology_file)
histology_file <- histology_file %>%
  dplyr::filter(short_histology %in% short_histology_of_interest)

# read gtf and filter to protein coding
cat('Reading Gencode file \n') 
gtf_file <- opt$gtf_file
gencode_gtf <- rtracklayer::import(con = gtf_file) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  dplyr::filter(!grepl("ENSG", gene_name), gene_type == "protein_coding") %>%
  unique()

# 1) read count data
cat('Filtering expression data \n')
count_file <- opt$count_file
count_mat <- readRDS(file = count_file)
# Convert matrix to data frame if it's a matrix
if(is.matrix(count_mat)) {
  count_mat <- as.data.frame(count_mat)
}
count_mat <- count_mat %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))

# filter expression count file to contain only protein coding genes
count_mat <- count_mat %>%
  dplyr::filter(rownames(count_mat) %in% gencode_gtf$gene_name)

# combine with histology to map sample id
count_mat <- reshape2::melt(as.matrix(count_mat),
                            varnames = c("Gene", "Kids_First_Biospecimen_ID"))
count_mat <- count_mat %>%
  dplyr::inner_join(histology_file %>% dplyr::select(Kids_First_Biospecimen_ID, sample_id),
                    by = 'Kids_First_Biospecimen_ID')

# get unique biospecimens per sample id
unique_ids <- count_mat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  dplyr::group_by(sample_id) %>%
  distinct(sample_id, .keep_all = T)
count_mat <- count_mat %>%
  dplyr::filter(
    sample_id %in% unique_ids$sample_id,
    Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID
  )
count_samples <- count_mat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix
count_mat <- count_mat %>%
  reshape2::acast(sample_id ~ Gene, value.var = "value")

# print dimensions
print(dim(count_mat))

# 2) Methylation
# read beta-values
cat('Reading beta values and subsetting \n')
methyl_file <- opt$methyl_file
methyl_data <- readRDS(file = file.path(methyl_file))
# Convert matrix to data frame if it's a matrix
if(is.matrix(methyl_data)) {
  methyl_data <- as.data.frame(methyl_data)
}
methyl_data <- methyl_data %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
methyl_data <- methyl_data[complete.cases(methyl_data), ]

# subset to unique ids
hist_methyl <- histology_file %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% colnames(methyl_data))
unique_ids <- hist_methyl %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::distinct(sample_id, .keep_all = T)
# Ensure methyl_data is a data frame before using select
if(is.matrix(methyl_data)) {
  methyl_data <- as.data.frame(methyl_data)
}
methyl_data <- methyl_data %>%
  dplyr::select(any_of(unique_ids$Kids_First_Biospecimen_ID))
hist_methyl <- hist_methyl %>%
  tibble::column_to_rownames("Kids_First_Biospecimen_ID")
hist_methyl <- hist_methyl[colnames(methyl_data), ]

# assign sample ids to columns to generate sample_id-probe matrix
stopifnot(identical(rownames(hist_methyl), colnames(methyl_data)))
colnames(methyl_data) <- hist_methyl$sample_id
methyl_data <- t(methyl_data) %>% as.data.frame()

# generate df of Kids_First_Biospecimen_ID and sample-id
methyl_samples <- hist_methyl %>%
  tibble::rownames_to_column("Kids_First_Biospecimen_ID") %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# print dimensions
print(dim(methyl_data))

# 3) Splice dataset
# read splice data
cat('Reading splice data and filtering \n')
splice_file <- opt$splice_file
splice_mat <- readRDS(splice_file)
# Convert matrix to data frame if it's a matrix
if(is.matrix(splice_mat)) {
  splice_mat <- as.data.frame(splice_mat)
}
splice_mat <- splice_mat %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
splice_mat <- reshape2::melt(as.matrix(splice_mat),
                             varnames = c("Splice_Variant", "Kids_First_Biospecimen_ID"))

# combine with histology to map sample id
splice_mat <- splice_mat %>%
  dplyr::inner_join(histology_file %>% dplyr::select(Kids_First_Biospecimen_ID, sample_id),
                    by = 'Kids_First_Biospecimen_ID')

# get unique biospecimens per sample id
unique_ids <- splice_mat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::distinct(sample_id, .keep_all = T)
splice_mat <- splice_mat %>%
  dplyr::filter(
    sample_id %in% unique_ids$sample_id,
    Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID
  )
splice_samples <- splice_mat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix
splice_mat <- splice_mat %>%
  acast(sample_id ~ Splice_Variant, value.var = "value")

# print dimensions
print(dim(splice_mat))

# subset to samples of interest (n = 152)
samples_of_interest <- intersect(rownames(count_mat), rownames(methyl_data))
samples_of_interest <- intersect(samples_of_interest, rownames(splice_mat))

# now final filter/transformation on samples of interest
cat('Performing feature selection \n')
# 1) RNA
# count_mat <- t(count_mat) %>% as.data.frame()
count_mat <- count_mat[samples_of_interest, ]
# remove genes with 0 counts across all samples
count_mat <- count_mat[, colSums(count_mat) > 0]

# Get top variable genes
num_genes = opt$num_features
keep <- apply(count_mat, 2, var)
keep <- keep[rev(order(keep))[1:num_genes]]
keep <- unique(c(names(keep)))
count_mat <- count_mat[, colnames(count_mat) %in% keep]

# rank transformation
count_mat <- datawizard::ranktransform(t(count_mat) %>% as.data.frame())
count_mat <- t(count_mat) %>% as.data.frame()
print(dim(count_mat)) # 1000
write_tsv(
  as.data.frame(count_mat) %>% rownames_to_column(),
  file = file.path(output_dir, "rna_data.tsv")
)

# 2) Methylation
methyl_data <- methyl_data[samples_of_interest, ]
methyl_data <- methyl_data[, colSums(is.na(methyl_data)) < nrow(methyl_data)]

# Get top variable methylation features
num_genes = opt$num_features
keep <- apply(methyl_data, 2, var)
keep <- keep[rev(order(keep))[1:num_genes]]
keep <- unique(c(names(keep)))
methyl_data <- methyl_data[, colnames(methyl_data) %in% keep]
print(dim(methyl_data)) # 1000
write_tsv(
  as.data.frame(methyl_data) %>% rownames_to_column(),
  file = file.path(output_dir, "methyl_data.tsv")
)

# 3) Splicing
# splice_mat <- t(splice_mat) %>% as.data.frame()
splice_mat <- splice_mat[samples_of_interest, ]
splice_mat <- splice_mat[, colSums(splice_mat) > 0]

# Get top variable splice features
num_features = opt$num_features
keep <- apply(splice_mat, 2, var)
keep <- keep[rev(order(keep))[1:num_features]]
keep <- unique(c(names(keep)))
splice_mat <- splice_mat[, colnames(splice_mat) %in% keep]
print(dim(splice_mat)) # 1000
write_tsv(
  as.data.frame(splice_mat) %>% rownames_to_column(),
  file = file.path(output_dir, "splice_data.tsv")
)

# final sample map
cat('Final step: Creating final sample map \n')
rna_samples <- count_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID")
methyl_samples <- methyl_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID")
splice_samples <- splice_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Splice" = "Kids_First_Biospecimen_ID")
sample_map <- rna_samples %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID_RNA) %>%
  inner_join(methyl_samples) %>%
  inner_join(splice_samples)
write_tsv(sample_map, file = file.path(output_dir, "samples_map.tsv"))
