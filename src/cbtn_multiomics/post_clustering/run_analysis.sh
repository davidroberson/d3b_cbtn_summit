# !/bin/bash
 
set -e
set -o pipefail
 
# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit
 
# Define directory and input files
data_dir="/sbgenomics/project-files/opc-v15"
data_prep_results="../data_preparation/results"
count_file="${data_prep_results}/rna_data.tsv"
methyl_file="${data_prep_results}/methyl_data.tsv"
splice_file="${data_prep_results}/splice_data.tsv"
cluster_file="../intNMF/results/intnmf_clusters.tsv"
output_dir="results/intnmf"
plots_dir="plots/intnmf"
histology_file="${data_dir}/v15/histologies.tsv"

# install packages
Rscript --vanilla 00-setup.R

# compare multi-modal derived clusters with molecular subtypes or subgroups
Rscript --vanilla 01-compare-classes.R \
--cluster_file $cluster_file \
--histology_file $histology_file \
--output_dir $output_dir

# feature and sample level heatmaps
Rscript --vanilla 02-heatmaps.R \
--cluster_file $cluster_file \
--histology_file $histology_file \
--rna_file $count_file \
--feature_scores_rna "../intNMF/results/feature_scores/feature_scores_rna.tsv" \
--methyl_file $methyl_file \
--feature_scores_methyl "../intNMF/results/feature_scores/feature_scores_methyl.tsv" \
--splice_file $splice_file \
--feature_scores_splice "../intNMF/results/feature_scores/feature_scores_splice.tsv" \
--plots_dir "${plots_dir}/heatmaps"

# generate balloon and corrplots
Rscript --vanilla 03-bubble-plots.R \
--cluster_file $cluster_file \
--histology_file $histology_file \
--plots_dir "${plots_dir}/bubble_plots"

# survival curves of Multi-modal clusters, RNA- and Methylation-derived subgroups (OS)
Rscript --vanilla 04-survival-curves_os.R \
--cluster_file $cluster_file \
--histology_file $histology_file \
--output_dir "${output_dir}/survival_os" \
--plots_dir "${plots_dir}/survival_os"

# survival curves of Multi-modal clusters, RNA- and Methylation-derived subgroups (EFS)
Rscript --vanilla 04-survival-curves_efs.R \
--cluster_file $cluster_file \
--histology_file $histology_file \
--plots_dir "${plots_dir}/survival_efs" \
--output_dir "${output_dir}/survival_efs"

# sankey plots showing relationship between Multi-modal clusters and other clinical variables
# Rscript --vanilla 05-sankey-plots.R \
# --cluster_file $cluster_file \
# --histology_file $histology_file \
# --plots_dir "${plots_dir}/sankey_plots"
