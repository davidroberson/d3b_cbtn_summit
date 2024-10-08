### Author: Komal S. Rathi, Adam Kraya

### Purpose

This module performs multi-modal clustering using the [IntNMF](https://cran.r-project.org/web/packages/IntNMF/index.html) R package on RNA, Methylation and Splicing data on the cohort of interest. 

### Run Analysis

```
# run full analysis
bash run_analysis.sh
```

***
`01-multi-modal-clustering-run.R`: Purpose of this script is to read input files generated in the `data_preparation` module and run IntNMF to identify the most optimal cluster fitting the dataset of interest.


#### Input

```
../data_preparation/results
├── methyl_data.tsv # methylation data used as input
├── rna_data.tsv # expression data used as input
├── splice_data.tsv # splice data used as input
└── samples_map.tsv # biospecimens + cohort identifiers for samples used for each modality 
```

#### Output

```
results/
├── feature_scores # IntNMF scores generated for each data modality
│ ├── feature_scores_methyl.tsv
│ ├── feature_scores_rna.tsv
│ └── feature_scores_splice.tsv
├── intnmf_best_fit.rds # output of nmf.mnnals for best fit (selected k)
├── intnmf_clusters.tsv # output clusters with per sample along with corresponding molecular subtype
└── intnmf_fit_all.rds # output of nmf.mnnals for all k values
```

#### Plots

```
plots
├── intnmf_consensus_plot.pdf # consensus plot of optimal k
└── intnmf_silhouette_plot.pdf # silhouette plot of optimal k
```
