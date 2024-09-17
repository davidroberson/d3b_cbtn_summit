#!/bin/bash
 
set -e
set -o pipefail
 
# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

echo "R_MAX_VSIZE=100Gb" > .Renviron

# update paths to new directory with subsetted data/matrices
data_dir="/sbgenomics/project-files/opc-v15"
histology_file="${data_dir}/v15/histologies.tsv"
short_histology="HGAT"
count_file="${data_dir}/v15/gene-counts-rsem-expected_count-collapsed.rds"
methyl_b_file="${data_dir}/v15/methyl-beta-values.rds"
functional_sites_splice_file="${data_dir}/v15/psi-se-primary.func.rds" # filtered to functional sites by Ammar Naqvi using v12
gtf_file="${data_dir}/v15/gencode.v39.primary_assembly.annotation.gtf.gz"

# prepare input files for multi-modal clustering
Rscript --vanilla 01-multi-modal-clustering-prepare-data.R \
--histology_file $histology_file \
--short_histology $short_histology \
--count_file $count_file \
--methyl_file $methyl_b_file \
--splice_file $functional_sites_splice_file \
--gtf_file $gtf_file \
--output_dir "results"
