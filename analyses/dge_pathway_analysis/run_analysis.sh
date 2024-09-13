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
gtf_file="${data_dir}/gencode.v39.primary_assembly.annotation.gtf.gz"
data_subset="../data_preparation/data"
count_file="${data_subset}/gene-counts-rsem-expected_count-collapsed.rds"
kegg_medicus_file="${data_dir}/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt"
cluster_file="../intNMF/results/intnmf_clusters.tsv"

# run DESeq2 analysis
Rscript --vanilla 01-dge_analysis_deseq.R \
--expr_mat $count_file \
--gtf_file $gtf_file \
--cluster_file $cluster_file \
--kegg_medicus_file $kegg_medicus_file \
--results_dir "results/intNMF/deseq" \
--plots_dir "plots/intNMF/deseq"
