### Author: Komal S. Rathi, Adam Kraya

### Purpose

Downstream analyses on methylation data identifying differentially methylated CpGs and regions and performing functional enrichment.

### Run Analysis

Full analysis can be run as follows:

```
bash run_analysis.sh
```

***
### Differentially methylated CpG site analysis

`01-limma_analysis.R`: The function of this script is to pull the full medulloblastoma methylation m-values dataset and filter the dataset into probes representing promoter region, gene-body (introns and exons) and promoter + gene-body combined. Next, differential methylation analyses is performed with the `limma` R package using a cluster-of-interest vs 'rest' approach (so 14 differential methylation analyses for the 14 clusters). 

Tests of interest (array.type = 'EPIC'):
a. gene body only
b. promoter only
c. gene body + promoter

#### Inputs

```
# intNMF derived clusters
../intNMF/results
└── intnmf_clusters.tsv

# methylation annotation file
../../data/v15
└── infinium.gencode.v39.probe.annotations.tsv.gz

# methylation m-values subsetted to cohort of interest
../../data/v15
└── methyl-m-values-hgat.rds
```

#### Outputs

TSV file of differentially methylation CPG sites (`FDR < 0.05`) for each cluster:

```
results
└── limma_output
    ├── gene_body_diffexpr_probes_per_cluster.tsv
    ├── genebody_promoter_diffexpr_probes_per_cluster.tsv
    └── promoter_diffexpr_probes_per_cluster.tsv
```

### Pathway enrichment of differentially methylated CpG sites

`02-dms_gsameth_analysis.R`: The purpose of this script is to use differentially methylated CpG sites obtained from limma analysis in `01-limma_analysis.R` and perform a pathway enrichment using `missMethyl::gsameth` function.

#### Inputs

```
# methylation m-values subsetted to cohort of interest
../../data/v15
└── methyl-m-values-hgat.rds

# differentially methylated CpG sites
results/limma_output
└── genebody_promoter_diffexpr_probes_per_cluster.tsv
```

#### Outputs

```
# table of pathway enrichment output (FDR < 0.1) on top 10000 differentially methylated CpG sites 
results
└── dms_gsameth_output
    └── hallmark
        └── genebody_promoter_gsameth_output_per_cluster.tsv

# plot of top 50 pathways (FDR < 0.1)
plots
└── dms_gsameth_output
    └── hallmark
        └── genebody_promoter_gsameth_pathways.pdf
```
