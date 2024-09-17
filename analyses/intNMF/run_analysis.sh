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
data_prep_results="../data_preparation/results"
samples_map="${data_prep_results}/samples_map.tsv"
count_file="${data_prep_results}/rna_data.tsv"
methyl_file="${data_prep_results}/methyl_data.tsv"
splice_file="${data_prep_results}/splice_data.tsv"

# run IntNMF multi-modal analysis
Rscript --vanilla 01-multi-modal-clustering-run.R \
--samples_map $samples_map \
--count_file $count_file \
--methyl_file $methyl_file \
--splice_file $splice_file \
--output_dir "results" \
--plots_dir "plots"
