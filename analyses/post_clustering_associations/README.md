### Author: Komal S. Rathi, Adam Kraya

### Purpose

This module generates survival plots for multi-modal derived clusters and balloon plots, corrplots, and sankey plots for associating multi-modal derived clusters to clinical variables like `molecular_subtype`, `CNS_region` etc. 

### Run Analysis

```
# run full analysis
bash run_analysis.sh
```

***
`01-compare-classes.R`: The purpose of this script is to perform chisq test of independence and adjusted rand index between multi-modal derived clusters and RNA-derived subtypes as well as with methylation-derived subclasses.

#### Inputs

```
# histology file
../../data/v15
└── histologies.tsv

# multi-modal derived clusters
../intNMF/results/
└── intnmf_clusters.tsv
```

#### Outputs

```
results
└── intnmf
    ├── ari_mm_vs_subtypes.txt # Adjusted rand index
    └── chisq_mm_vs_subtypes.txt # Chisq test of independence
```

***
`02-heatmaps.R`: The purpose of this script is to generate feature and sample level heatmaps for multi-modal derived clusters.

#### Inputs

```
# histology file
../../data/v15
└── histologies.tsv

# filtered and transformed files used for multi-modal clustering
../data_preparation/results
├── methyl_data.tsv
├── rna_data.tsv
└── splice_data.tsv

# feature-level intNMF weights for each data modality
../intNMF/results/feature_scores
├── feature_scores_methyl.tsv
├── feature_scores_rna.tsv
└── feature_scores_splice.tsv

# multi-modal derived clusters
../intNMF/results
└── intnmf_clusters.tsv
```

#### Outputs

```
plots/intnmf/heatmaps
├── feature_level_heatmaps.pdf
└── sample_level_heatmaps.pdf
```

***

`03-bubble-plots.R`: The purpose of this script is to generate bubble plots like balloon plots and corrplots to identify the association between multi-modal derived clusters and  RNA-derived molecular subtypes and methylation-derived subclasses.

#### Inputs

```
# histology file
../../data/v15
└── histologies.tsv

# multi-modal derived clusters
../intNMF/results
└── intnmf_clusters.tsv
```

#### Outputs

```
plots/intnmf/bubble_plots
├── mm_clusters_vs_dkfz_v11_methylation_subclass_balloonplot.pdf
├── mm_clusters_vs_dkfz_v11_methylation_subclass_corrplot.pdf
├── mm_clusters_vs_dkfz_v12_methylation_subclass_balloonplot.pdf
├── mm_clusters_vs_dkfz_v12_methylation_subclass_corrplot.pdf
├── mm_clusters_vs_molsubtype_balloonplot.pdf
└── mm_clusters_vs_molsubtype_corrplot.pdf
```
***
`04-survival-curves_os.R`: Purpose of this script is to generate overall survival curves of multi-modal derived clusters, RNA-derived molecular subtypes and methylation-derived subclasses.

#### Inputs

```
# histology file
../../data/v15
└── histologies.tsv

# multi-modal derived clusters
../intNMF/results
└── intnmf_clusters.tsv
```

#### Outputs

```
# risk scores and summaries
results/intnmf/survival_os
├── coxph_risk_score_OS.txt
└── coxph_summary_OS.txt

# overall survival plots
plots/intnmf/survival_os
├── coxph_summary_OS.pdf
├── survival_methyl_subtype_v11.pdf
├── survival_methyl_subtype_v12.pdf
├── survival_mm_clusters.pdf
└── survival_molsubtype.pdf
```
***
`04-survival-curves_efs.R`: Purpose of this script is to generate event-free survival curves of multi-modal derived clusters, RNA-derived molecular subtypes and methylation-derived subclasses.

#### Inputs

```
# histology file
../../data/v15
└── histologies.tsv

# multi-modal derived clusters
../intNMF/results
└── intnmf_clusters.tsv
```

#### Outputs

```
# risk scores and summaries
results/intnmf/survival_efs
├── coxph_risk_score_EFS.txt
└── coxph_summary_EFS.txt

# event-free survival plots
plots/intnmf/survival_efs
├── coxph_summary_EFS.pdf
├── survival_methyl_subtype_v11.pdf
├── survival_methyl_subtype_v12.pdf
├── survival_mm_clusters.pdf
└── survival_molsubtype.pdf
```

***
`05-sankey-plots.R`:  The purpose of this script is to generate sankey plots showing distribution between:

1. RNA-derived molecular subtypes, multi-modal derived clusters and CNS region
2. Methylation-derived subclass, multi-modal derived clusters and CNS region
3. RNA-derived molecular subtypes, multi-modal derived clusters and EFS_event_type
4. Methylation-derived subclass, multi-modal derived clusters and EFS_event_type

#### Inputs

```
# for sankey plots, we need EFS_event_type which is present in v15 histology file
../../data/v15
└── histologies.tsv

# multi-modal derived clusters
../intNMF/results
└── intnmf_clusters.tsv
```

#### Output

```
plots/intnmf/sankey_plots
├── dkfz_v12_methylation_subclass_vs_mm_clusters_cns_region_sankey.pdf
├── dkfz_v12_methylation_subclass_vs_mm_clusters_efs_event_sankey.pdf
├── molecular_subtype_vs_mm_clusters_cns_region_sankey.pdf
└── molecular_subtype_vs_mm_clusters_efs_event_sankey.pdf
```
