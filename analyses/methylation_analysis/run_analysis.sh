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

# Define directory and input files
data_dir="/sbgenomics/project-files/opc-v15"
count_file="${data_dir}/v15/gene-counts-rsem-expected_count-collapsed.rds"
methyl_m_file="${data_dir}/v15/methyl-m-values.rds" 
cluster_file="../intNMF/results/intnmf_clusters.tsv"
methyl_annot_file="${data_dir}/v15/infinium.gencode.v39.probe.annotations.tsv.gz"

# 1) limma analysis to identify differentially expressed CpG sites
Rscript --vanilla 01-limma_analysis.R \
--methyl_mat $methyl_m_file \
--methyl_annot $methyl_annot_file \
--cluster_file $cluster_file \
--output_dir "results/limma_output"

# 2) gsameth pathway enrichment on differentially expressed CpG sites from (genebody + promoter) obtained from 01-limma_analysis.R
# using REACTOME
Rscript --vanilla 02-dms_gsameth_analysis.R \
--methyl_mat $methyl_m_file \
--diffexpr_sites "results/limma_output/genebody_promoter_diffexpr_probes_per_cluster.tsv" \
--msigdb "reactome" \
--prefix "genebody_promoter" \
--output_dir "results/dms_gsameth_output/reactome" \
--plots_dir "plots/dms_gsameth_output/reactome"

# using HALLMARK
Rscript --vanilla 02-dms_gsameth_analysis.R \
--methyl_mat $methyl_m_file \
--diffexpr_sites "results/limma_output/genebody_promoter_diffexpr_probes_per_cluster.tsv" \
--msigdb "hallmark" \
--prefix "genebody_promoter" \
--output_dir "results/dms_gsameth_output/hallmark" \
--plots_dir "plots/dms_gsameth_output/hallmark"

# 2) methreg analysis on differentially methylated CpG sites from (genebody + promoter) obtained from 01-limma_analysis.R
Rscript --vanilla 02-methreg_analysis.R \
--methyl_mat $methyl_m_file \
--cluster_file $cluster_file \
--expr_mat $count_file \
--diffexpr_sites "results/limma_output/genebody_promoter_diffexpr_probes_per_cluster.tsv" \
--output_dir "results/methreg_output" \
--plots_dir "plots/methreg_output"

# 3) identify differentially methylation regions using DMRcate::dmrcate and pathway enrichment using missMethyl::gsaregion
# using REACTOME
Rscript --vanilla 03-dmr_gsaregion_analysis.R \
--methyl_mat $methyl_m_file \
--methyl_annot $methyl_annot_file \
--cluster_file $cluster_file \
--msigdb "reactome" \
--output_dir "results/dmr_gsaregion_output/reactome" \
--plots_dir "plots/dmr_gsaregion_output/reactome"

# using HALLMARK
Rscript --vanilla 03-dmr_gsaregion_analysis.R \
--methyl_mat $methyl_m_file \
--methyl_annot $methyl_annot_file \
--cluster_file $cluster_file \
--msigdb "hallmark" \
--output_dir "results/dmr_gsaregion_output/hallmark" \
--plots_dir "plots/dmr_gsaregion_output/hallmark"


