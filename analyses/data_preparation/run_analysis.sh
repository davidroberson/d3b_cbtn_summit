#!/bin/bash
 
set -e
set -o pipefail
 
# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# update paths to new directory with subsetted data/matrices
data_dir="../../data"
histology_file="${data_dir}/v15/histologies.tsv"
short_histology="HGAT"
count_file="${data_dir}/v15/gene-counts-rsem-expected_count-collapsed.rds"
cnv_gainloss_file="${data_dir}/v15/All.gainloss.txt.gz"
snv_file="${data_dir}/v15/snv-consensus-plus-hotspots.maf.tsv.gz"
methyl_b_file="${data_dir}/v15/methyl-beta-values.rds"
functional_sites_splice_file="${data_dir}/v15/splice_events_pan_cancer_functional_filter.rds" # filtered to functional sites by Ammar Naqvi using v12

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
