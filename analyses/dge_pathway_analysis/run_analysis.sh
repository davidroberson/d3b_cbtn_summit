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
gtf_file="${data_dir}/v15/gencode.v39.primary_assembly.annotation.gtf.gz"
count_file="${data_dir}/v15/gene-counts-rsem-expected_count-collapsed.rds"
cluster_file="../intNMF/results/intnmf_clusters.tsv"

# run DESeq2 analysis
Rscript --vanilla 01-dge_analysis_deseq.R \
--expr_mat $count_file \
--gtf_file $gtf_file \
--cluster_file $cluster_file \
--results_dir "results/intNMF/deseq" \
--plots_dir "plots/intNMF/deseq"
