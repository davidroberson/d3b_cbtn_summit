#!/bin/bash
 
set -e
set -o pipefail
 
# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit
 
# Define directory and input files
data_dir="../../data"
histology_file="${data_dir}/v15/histologies.tsv"
short_histology="HGAT"
count_file="${data_dir}/v15/gene-counts-rsem-expected_count-collapsed.rds"
tpm_file="${data_dir}/v15/gene-expression-rsem-tpm-collapsed.rds"
snv_file="${data_dir}/v15/snv-consensus-plus-hotspots.maf.tsv.gz"
methyl_m_file="${data_dir}/v15/methyl-m-values.rds" 
methyl_b_file="${data_dir}/v15/methyl-beta-values.rds" 
cnv_gainloss_file="${data_dir}/20230309_release.All.gainloss.txt.gz" # from data assembly 20230309 release 
rmats_splice_file="${data_dir}/splice-events-rmats.tsv.gz" # soft-linked to pbta-splicing/data/v5/splice-events-rmats.tsv.gz
functional_sites_splice_file="${data_dir}/splice_events_pan_cancer_functional_filter.rds" # filtered to functional sites by Ammar Naqvi using v12
gtf_file="${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz"
cancer_genes="${data_dir}/cancer_gene_list.rds"

# subset data using short histology for multi-modal clustering
Rscript --vanilla 01-subset-data.R \
--histology_file $histology_file \
--short_histology $short_histology \
--count_file $count_file \
--tpm_file $tpm_file \
--snv_file $snv_file \
--methyl_m_file $methyl_m_file \
--methyl_b_file $methyl_b_file \
--splice_file $functional_sites_splice_file \
--cnv_file $cnv_gainloss_file \
--output_dir "data"

# update paths to new directory with subsetted data/matrices
data_dir="data"
histology_file="${data_dir}/histologies.tsv"
short_histology="HGAT"
count_file="${data_dir}/gene-counts-rsem-expected_count-collapsed.rds"
cnv_gainloss_file="${data_dir}/All.gainloss.txt.gz"
snv_file="${data_dir}/snv-consensus-plus-hotspots.maf.tsv.gz"
methyl_b_file="${data_dir}/methyl-beta-values.rds"
functional_sites_splice_file="${data_dir}/splice_events_pan_cancer_functional_filter.rds" # filtered to functional sites by Ammar Naqvi using v12

# prepare input files for multi-modal clustering
Rscript --vanilla 03-multi-modal-clustering-prepare-data.R \
--histology_file $histology_file \
--short_histology $short_histology \
--cancer_genes $cancer_genes \
--count_file $count_file \
--cnv_gainloss_file $cnv_gainloss_file \
--snv_file $snv_file \
--methyl_file $methyl_b_file \
--splice_file $functional_sites_splice_file \
--gtf_file $gtf_file \
--output_dir "results"
