#!/usr/bin/env Rscript

# Script to create minimal test data for the multi-modal clustering workflow
# Simple version without tidyverse dependency

# Create directory structure
dir.create("test_data", showWarnings = FALSE, recursive = TRUE)

# ===== 1. Create minimal histology file =====
sample_ids <- paste0("sample", 1:20)
bs_ids <- paste0("BS", 1:20)
histology <- data.frame(
  sample_id = sample_ids,
  Kids_First_Biospecimen_ID = bs_ids,
  short_histology = rep("HGAT", 20),
  molecular_subtype = sample(c("Subtype1", "Subtype2", "Subtype3"), 20, replace = TRUE),
  dkfz_methylation_class = sample(c("Class1", "Class2", "Class3"), 20, replace = TRUE),
  sex = sample(c("Male", "Female"), 20, replace = TRUE),
  age_at_diagnosis = runif(20, 0, 18),
  OS_days = runif(20, 100, 2000),
  OS_status = sample(c("Alive", "Dead"), 20, replace = TRUE),
  EFS_days = runif(20, 100, 2000),
  EFS_status = sample(c("Event", "Censored"), 20, replace = TRUE),
  CNS_region = sample(c("Region1", "Region2", "Region3"), 20, replace = TRUE)
)

# Write histology file
write.table(histology, "test_data/histologies.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# ===== 2. Create minimal RNA-seq data =====
# Create 1000 genes x 20 samples matrix
gene_names <- paste0("gene", 1:1000)
rna_data <- matrix(rpois(1000 * 20, lambda = 100), nrow = 1000)
rownames(rna_data) <- gene_names
colnames(rna_data) <- histology$Kids_First_Biospecimen_ID

# Save as RDS
saveRDS(rna_data, "test_data/gene-counts-rsem-expected_count-collapsed.rds")

# ===== 3. Create minimal methylation data =====
# Create 1000 CpGs x 20 samples matrix
cpg_names <- paste0("cg", 1:1000)
methyl_data <- matrix(runif(1000 * 20, 0, 1), nrow = 1000)
rownames(methyl_data) <- cpg_names
colnames(methyl_data) <- histology$Kids_First_Biospecimen_ID

# Save as RDS
saveRDS(methyl_data, "test_data/methyl-beta-values.rds")

# ===== 4. Create minimal splicing data =====
# Create 1000 junctions x 20 samples matrix
junction_names <- paste0("junction", 1:1000)
splice_data <- matrix(runif(1000 * 20, 0, 1), nrow = 1000)
rownames(splice_data) <- junction_names
colnames(splice_data) <- histology$Kids_First_Biospecimen_ID

# Save as RDS
saveRDS(splice_data, "test_data/psi-se-primary.func.rds")

# ===== 5. Create minimal GTF file =====
gtf_content <- paste0(
  "chr1\thavana\tgene\t1\t1000\t.\t+\t.\tgene_id \"gene", 1:1000, "\"; gene_name \"gene", 1:1000, "\"; gene_type \"protein_coding\";\n", 
  collapse = ""
)

# Write GTF content to a temporary file
gtf_file <- "test_data/gencode.v39.primary_assembly.annotation.gtf"
writeLines(gtf_content, gtf_file)

# Compress the GTF file
system(paste("gzip", gtf_file))

# ===== 6. Create input YAML file =====
yaml_content <- paste0(
'# Test input parameters for multi-omic clustering workflow
histology_file:
  class: File
  path: test_data/histologies.tsv
short_histology: "HGAT"
count_file:
  class: File
  path: test_data/gene-counts-rsem-expected_count-collapsed.rds
methyl_file:
  class: File
  path: test_data/methyl-beta-values.rds
splice_file:
  class: File
  path: test_data/psi-se-primary.func.rds
gtf_file:
  class: File
  path: test_data/gencode.v39.primary_assembly.annotation.gtf.gz
num_features: 100
max_k: 3
'
)

# Write YAML file
writeLines(yaml_content, "test_data/test_inputs.yaml")

cat("Test data created successfully!\n")
cat("To test the workflow:\n")
cat("cwltool multi_modal_clustering_workflow.cwl test_data/test_inputs.yaml\n")