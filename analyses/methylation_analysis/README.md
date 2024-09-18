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
../intNMF/results/intnmf_clusters.tsv

# methylation annotation file
../../data/v15/infinium.gencode.v39.probe.annotations.tsv.gz

# methylation m-values subsetted to cohort of interest
../data_preparation/data
└── methyl-m-values.rds
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
../data_preparation/data
└── methyl-m-values.rds

# differentially methylated CpG sites
results/limma_output
└── genebody_promoter_diffexpr_probes_per_cluster.tsv
```

#### Outputs

```
# table of pathway enrichment output (FDR < 0.1) on top 10000 differentially methylated CpG sites 
results
└── dms_gsameth_output
    ├── hallmark
    │   └── genebody_promoter_gsameth_output_per_cluster.tsv
    └── reactome
        └── genebody_promoter_gsameth_output_per_cluster.tsv

# plot of top 50 pathways (FDR < 0.1)
plots
└── dms_gsameth_output
    ├── hallmark
    │   └── genebody_promoter_gsameth_pathways.pdf
    └── reactome
        └── genebody_promoter_gsameth_pathways.pdf
```

### Differential methylated region analysis (DMRcate) + Pathway Enrichment (gsaregion)

`03-dmr_gsaregion_analysis.R`:  The function of this script is to pull the full medulloblastoma methylation m-values dataset and filter the dataset into probes representing promoter region, gene-body (introns and exons) and promoter + gene-body combined. Next, differential region-level methylation analyses is performed with `DMRcate::dmrcate` using a cluster-of-interest vs 'rest' approach, and specifying 'EPIC' array. Finally, `missMethyl::gsaregion` analysis (setting sig.genes = TRUE) is performed on the the resulting differentially methylated regions to determine cluster-specific pathway differences. 

Tests of interest (array.type = 'EPIC'):
a. gene body only
b. promoter only
c. gene body + promoter

#### Inputs

```
# intNMF derived clusters
../intNMF/results/intnmf_clusters.tsv

# methylation annotation file
../../data/v15/infinium.gencode.v39.probe.annotations.tsv.gz

# methylation m-values subsetted to cohort of interest
../data_preparation/data
└── methyl-m-values.rds
```

TSV file of all significant pathways (`FDR < 0.05`) for each cluster:

```
results
└── dmr_gsaregion_output
    ├── hallmark
    │   ├── gene_body_gsaregion_output_per_cluster.tsv
    │   ├── genebody_promoter_gsaregion_output_per_cluster.tsv
    │   └── promoter_gsaregion_output_per_cluster.tsv
    └── reactome
        ├── gene_body_gsaregion_output_per_cluster.tsv
        ├── genebody_promoter_gsaregion_output_per_cluster.tsv
        └── promoter_gsaregion_output_per_cluster.tsv
```

Barplots of top 50 pathways (`FDR < 0.05`) enriched in each cluster:

```
plots
└── dmr_gsaregion_output
    ├── hallmark
    │   ├── gene_body_gsaregion_pathways.pdf
    │   ├── genebody_promoter_gsaregion_pathways.pdf
    │   └── promoter_gsaregion_pathways.pdf
    └── reactome
        ├── gene_body_gsaregion_pathways.pdf
        ├── genebody_promoter_gsaregion_pathways.pdf
        └── promoter_gsaregion_pathways.pdf
```

