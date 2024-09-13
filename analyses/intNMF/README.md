### Author: Komal S. Rathi, Adam Kraya

### Purpose

This module performs multi-modal clustering using the [IntNMF](https://cran.r-project.org/web/packages/IntNMF/index.html) R package on RNA, CNV, SNV, Methylation and Splicing data on the cohort of interest. 

### Run Analysis

```
# run full analysis
bash run_analysis.sh
```

***
`01-multi-modal-clustering-run.R`: Purpose of this script is to read input files generated in the `data_preparation` module and run IntNMF to identify the most optimal cluster fitting the dataset of interest.

#### Selection of optimal cluster

From this comment: https://github.com/d3b-center/bixu-tracker/issues/1704#issuecomment-1480304368
For each attempted k, the cluster values were mapped to the samples for each mode of data (e.g. CNV, SNV, Methylation, RNA and Splicing), the fpc stats for each mode were computed, sum of the `average.between` and `average.within` across the modes of data at each k were taken, then the difference of those two summed values at each k was taken to select the optimal cluster number (e.g. the k with the largest difference between the two).

Using this method, the samples were classified into `14 clusters`.

#### Input

```
../data_preparation/results
├── cnv_data.tsv # cnv data used as input
├── methyl_data.tsv # methylation data used as input
├── norm_counts.tsv # expression data used as input
├── snv_data.tsv # snv data used as input
├── splice_data.tsv # splice data used as input
└── samples_map.tsv # biospecimens + cohort identifiers for samples used for each modality 
```

#### Output

```
results
├── intnmf_fit_all.rds # output of nmf.mnnals for all k values
├── intnmf_clusterstats.tsv # cluster stats across all k-values for each modality
├── intnmf_best_fit.rds # output of nmf.mnnals for best fit (selected k)
├── feature_scores # IntNMF scores generated for each data modality
│   ├── feature_scores_cnv.tsv
│   ├── feature_scores_methyl.tsv
│   ├── feature_scores_rna.tsv
│   ├── feature_scores_snv.tsv
│   └── feature_scores_splice.tsv
└── intnmf_clusters.tsv # output clusters with per sample along with corresponding molecular subtype

plots
├── intnmf_consensus_plot.pdf # consensus plot of optimal k
└── intnmf_silhouette_plot.pdf # silhouette plot of optimal k
```
