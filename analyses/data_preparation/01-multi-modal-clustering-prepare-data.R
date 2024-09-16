# prepare files for multi-modal clustering
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(datawizard)
  library(reshape2)
})

# parse command line options
option_list <- list(
  make_option(c("--histology_file"), type = "character", help = "Histology file (.tsv)"),
  make_option(c("--short_histology"), type = "character", help = "Short histology of interest"),
  make_option(c("--cancer_genes"), type = "character", help = "cancer gene list (.rds)"),
  make_option(c("--count_file"), type = "character", help = "Expression matrix, preferably counts (.rds)"),
  make_option(c("--cnv_gainloss_file"), type = "character", help = "CNVkit gainloss file (.tsv)"),
  make_option(c("--snv_file"), type = "character", help = "SNV consensus (.maf | tsv.gz)"),
  make_option(c("--methyl_file"), type = "character", help = "Methylation values, preferably beta values (.rds)"),
  make_option(c("--splice_file"), type = "character", help = "PSI values (.rds)"),
  make_option(c("--gtf_file"), type = "character", help = "Gencode gtf file"),
  make_option(c("--output_dir"), type = "character", help = "path to results directory")
)
opt <- parse_args(OptionParser(option_list = option_list, add_help_option = TRUE))
short_histology_of_interest <- opt$short_histology

# output directory
output_dir <- opt$output_dir
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path("utils", "filter_cnv.R"))

# cancer genes
cancer_genes <- readRDS(opt$cancer_genes)

# read histology file and filter to short histology of interest
histology_file <- opt$histology_file
histology_file <- read_tsv(file = histology_file)
histology_file <- histology_file %>%
  dplyr::filter(short_histology %in% short_histology_of_interest)

# read gtf and filter to protein coding
gtf_file <- opt$gtf_file
gencode_gtf <- rtracklayer::import(con = gtf_file) %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  dplyr::filter(!grepl("ENSG", gene_name), gene_type == "protein_coding") %>%
  unique()

# 1) read count data
count_file <- opt$count_file
count_mat <- readRDS(file = count_file)
count_mat <- count_mat %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))

# filter expression count file to contain only protein coding genes
count_mat <- count_mat %>%
  dplyr::filter(rownames(count_mat) %in% gencode_gtf$gene_name)

# combine with histology to map sample id
count_mat <- melt(as.matrix(count_mat),
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

# 2) CNV
cnv_gainloss_file <- opt$cnv_gainloss_file
cnv_dat <- data.table::fread(file = cnv_gainloss_file)
cnv_dat <- cnv_dat %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "BS_ID") %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% histology_file$Kids_First_Biospecimen_ID) %>%
  dplyr::select(Kids_First_Biospecimen_ID, gene, log2) %>%
  dplyr::inner_join(
    histology_file %>% dplyr::select(
      Kids_First_Biospecimen_ID,
      sample_id,
      tumor_ploidy,
      tumor_fraction
    ),
    by = "Kids_First_Biospecimen_ID"
  ) %>%
  unique()

# get unique biospecimens per sample id
unique_ids <- cnv_dat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::distinct(sample_id, .keep_all = T)
cnv_dat <- cnv_dat %>%
  dplyr::filter(
    sample_id %in% unique_ids$sample_id,
    Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID
  )
cnv_samples <- cnv_dat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix
# function to adjust copy number and status
adjust_cn <- function(x) {
  # get tumor fraction and ploidy
  tumor_fraction <- unique(na.omit(x$tumor_fraction))
  tumor_ploidy <- unique(na.omit(x$tumor_ploidy))
  
  if (length(tumor_fraction) == 1 & length(tumor_ploidy) == 1) {
    # calculate adjusted copy number if tumor fraction and ploidy info is available
    x$adjusted_cn <- (((2 ^ (x$log2) - (
      1 - tumor_fraction
    )) * tumor_ploidy) / tumor_fraction) - 0.5
    x$adjusted_cn <- round(x$adjusted_cn)
    x$adjusted_status <- ifelse(x$adjusted_cn == 0,
                                "Complete Loss",
                                ifelse(
                                  x$adjusted_cn == 1,
                                  "Loss",
                                  ifelse(
                                    x$adjusted_cn %in% c(tumor_ploidy + 1:9),
                                    "Gain",
                                    ifelse(x$adjusted_cn >= 10, "Amplification", "Neutral")
                                  )
                                ))
    
    # replace old columns with new ones
    x <- x %>%
      dplyr::rename("status" = "adjusted_status", # rename new columns
                    "copy_number" = "adjusted_cn")
    
  }
}

# apply function to all samples in the consensus file
cnv_dat <- plyr::ddply(
  .data = cnv_dat,
  .variables = "sample_id",
  .fun = function(x)
    adjust_cn(x = x)
)
cnv_dat <- cnv_dat %>%
  dplyr::filter(status != "Neutral") %>%
  dplyr::select(-c(copy_number)) %>%
  unique()

# subset to cancer genes
cnv_dat <- cnv_dat %>%
  dplyr::rename("hgnc_symbol" = "gene") %>%
  filter_cnv(myCancerGenes = cancer_genes)
cnv_dat <- acast(
  cnv_dat,
  sample_id ~ hgnc_symbol,
  value.var = "log2",
  fill = 0,
  fun.aggregate = max
)

cnv_dat <- abs(cnv_dat)

# print dimensions
print(dim(cnv_dat))

# 3) SNV
snv_file <- opt$snv_file
snv_dat <- data.table::fread(file = snv_file)
maf_nonsynonymous <- c(
  "Missense_Mutation",
  "Frame_Shift_Del",
  "In_Frame_Ins",
  "Frame_Shift_Ins",
  "Splice_Site",
  "Nonsense_Mutation",
  "In_Frame_Del",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)
snv_dat <- snv_dat %>%
  dplyr::filter(
    Variant_Classification %in% maf_nonsynonymous,
    Tumor_Sample_Barcode %in% histology_file$Kids_First_Biospecimen_ID
  )

# combine with histology to map sample id
snv_dat <- snv_dat %>%
  dplyr::rename("Kids_First_Biospecimen_ID" = "Tumor_Sample_Barcode") %>%
  dplyr::inner_join(histology_file %>% dplyr::select(Kids_First_Biospecimen_ID, sample_id),
             by = 'Kids_First_Biospecimen_ID')

# get unique biospecimens per sample id
unique_ids <- snv_dat %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID) %>%
  unique() %>%
  dplyr::arrange(Kids_First_Biospecimen_ID) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::distinct(sample_id, .keep_all = T)
snv_dat <- snv_dat %>%
  dplyr::filter(
    sample_id %in% unique_ids$sample_id,
    Kids_First_Biospecimen_ID %in% unique_ids$Kids_First_Biospecimen_ID
  )
snv_samples <- snv_dat %>%
  dplyr::select(Kids_First_Biospecimen_ID, sample_id) %>%
  unique()

# convert to matrix
snv_dat <- acast(snv_dat, sample_id ~ Hugo_Symbol, value.var = "Variant_Classification")
snv_dat[snv_dat > 0] <- 1 # anything > 0 should be coded as 1

# print dimensions
print(dim(snv_dat))

# 4) Methylation
# read beta-values
methyl_file <- opt$methyl_file
methyl_data <- readRDS(file = file.path(methyl_file))
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

# 5) Splice dataset
# read splice data
splice_file <- opt$splice_file
splice_mat <- readRDS(splice_file)
splice_mat <- splice_mat %>%
  dplyr::select(any_of(histology_file$Kids_First_Biospecimen_ID))
splice_mat <- melt(as.matrix(splice_mat),
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
samples_of_interest <- intersect(rownames(count_mat), rownames(cnv_dat))
samples_of_interest <- intersect(samples_of_interest, rownames(snv_dat))
samples_of_interest <- intersect(samples_of_interest, rownames(methyl_data))
samples_of_interest <- intersect(samples_of_interest, rownames(splice_mat))

# now final filter/transformation on samples of interest

# 1) RNA
# count_mat <- t(count_mat) %>% as.data.frame()
count_mat <- count_mat[samples_of_interest, ]
# remove genes with 0 counts across all samples
count_mat <- count_mat[, colSums(count_mat) > 0]

# top 1000 most variable genes
num_genes = 1000
keep <- apply(count_mat, 2, var)
keep <- keep[rev(order(keep))[1:num_genes]]
keep <- unique(c(names(keep)))
count_mat <- count_mat[, colnames(count_mat) %in% keep]

# rank transformation
count_mat <- ranktransform(t(count_mat) %>% as.data.frame())
count_mat <- t(count_mat) %>% as.data.frame()
print(dim(count_mat)) # 1000
write_tsv(
  as.data.frame(count_mat) %>% rownames_to_column(),
  file = file.path(output_dir, "rna_data.tsv")
)

# 2) CNV
# filter for standard deviation (too many features)
cnv_dat <- cnv_dat[samples_of_interest, ]
cnv_dat <- as.data.frame(cnv_dat)
print(dim(cnv_dat)) # 1193
write_tsv(as.data.frame(cnv_dat) %>% rownames_to_column(),
          file = file.path(output_dir, "cnv_data.tsv"))

# 3) SNV
# filter out genes with low mutation rate
snv_dat <- snv_dat[samples_of_interest, ]
mut.rate <- apply(X = snv_dat, MARGIN = 2, FUN = mean)
snv_dat <- snv_dat[, which(mut.rate > 0.02)]
print(dim(snv_dat)) # 170
write_tsv(as.data.frame(snv_dat) %>% rownames_to_column(),
          file = file.path(output_dir, "snv_data.tsv"))

# 4) Methylation
methyl_data <- methyl_data[samples_of_interest, ]
methyl_data <- methyl_data[, colSums(is.na(methyl_data)) < nrow(methyl_data)]

# top 1000 most variable features
num_genes = 1000
keep <- apply(methyl_data, 2, var)
keep <- keep[rev(order(keep))[1:num_genes]]
keep <- unique(c(names(keep)))
methyl_data <- methyl_data[, colnames(methyl_data) %in% keep]
print(dim(methyl_data)) # 1000
write_tsv(
  as.data.frame(methyl_data) %>% rownames_to_column(),
  file = file.path(output_dir, "methyl_data.tsv")
)

# 5) Splicing
# splice_mat <- t(splice_mat) %>% as.data.frame()
splice_mat <- splice_mat[samples_of_interest, ]
splice_mat <- splice_mat[, colSums(splice_mat) > 0]

# top 1000 most variable features
num_features = 1000
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
rna_samples <- count_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_RNA" = "Kids_First_Biospecimen_ID")
methyl_samples <- methyl_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Methyl" = "Kids_First_Biospecimen_ID")
cnv_samples <- cnv_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_CNV" = "Kids_First_Biospecimen_ID")
snv_samples <- snv_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_SNV" = "Kids_First_Biospecimen_ID")
splice_samples <- splice_samples %>%
  dplyr::filter(sample_id %in% samples_of_interest) %>%
  dplyr::rename("Kids_First_Biospecimen_ID_Splice" = "Kids_First_Biospecimen_ID")
sample_map <- rna_samples %>%
  dplyr::select(sample_id, Kids_First_Biospecimen_ID_RNA) %>%
  inner_join(cnv_samples) %>%
  inner_join(snv_samples) %>%
  inner_join(methyl_samples) %>%
  inner_join(splice_samples)
write_tsv(sample_map, file = file.path(output_dir, "samples_map.tsv"))
