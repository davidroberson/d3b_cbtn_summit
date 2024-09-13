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
../../data/v12/infinium.gencode.v39.probe.annotations.tsv.gz

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

### Functional prioritization of differentially methylated CpG sites

`02-methreg_analysis.R`: The R package `MethReg` uses triplet information from `DNA methylation + RNA expression data + transcription factor (TF) binding sites` for functional prioritization of significantly differentially methylated sites. 

The input for this script is *matched* `methylation m-values` and `RNA expression matrix (raw counts)` matrix. 

- There are two types of workflows: supervised and unsupervised.

  1. The `unsupervised workflow` does not account for any design information (i.e. cannot differentiate groups in the matrix) so the solution would be running on each cluster individually by subsetting to samples of interest and running the workflow.  
  2. For the `supervised workflow`, we can use the cluster-specific DMS or DMRs as input. 

**Note**: We chose to use the DMSites from the `gene body + promoter` output from step 1) i.e. the output of `01-limma_analysis.R` to run the supervised workflow. 

- Creating the triplets (DNA methylation + RNA expression + TF binding sites) can be done in three ways as described  [here](https://www.bioconductor.org/packages/release/bioc/vignettes/MethReg/inst/doc/MethReg.html#methreg-workflow). 

    1.  Mapping the region to the closest gene (`target.method = "genes.promoter.overlap"`).
    2.  Mapping the region to a specific number of genes upstream down/upstream of the region (`target.method = "nearby.genes"`) 
    3.  Mapping the region to all the genes within a window (default size = 500 kbp around the region, i.e. +- 250 kbp from start or end of the region) (`target.method = "window"`)  
 
Linear model is applied where methylation levels are represented as binary values for high and low values.

**Note**: We used `Mapping the region to the closest gene` because it takes the least amount of runtime. The other two methods show a runtime of 25-30h.

#### Inputs

```
# methylation m-values subsetted to cohort of interest
../data_preparation/data
└── methyl-m-values.rds

# gene expression 
../data_preparation/data
└── gene-counts-rsem-expected_count-collapsed.rds

# intNMF derived clusters
../intNMF/results/intnmf_clusters.tsv

# differentially methylated CpG sites
results/limma_output
└── genebody_promoter_diffexpr_probes_per_cluster.tsv
```

#### Outputs

- Output of the workflow:

  1.  Table of significant triplets with the role of TF in the target gene expression (repressor or activator) and the role of DNA methylation on TF (enhancing or attenuating). 
  2. Table of fitting the linear model with p-value and estimate of 1) direct effect of DNAm 2) direct effect of TF 3) synergistic effect of DNAm + TF on gene expression.
  3. Associated scatter plots of top 5 triplets per cluster. 

Here is the output of supervised workflow with DMSites obtained after running `01-limma_analysis.R` + triplet creation using mapping TF binding sites to the closest gene:

1. `triplet_nearest_gene_interactions.tsv`: Table of fitting the linear model with p-value and estimate of 1) direct effect of DNAm 2) direct effect of TF 3) synergistic effect of DNAm + TF on gene expression.

2. `triplet_nearest_gene_interactions_stratified_model.tsv`: Table of significant triplets with the role of TF in the target gene expression (repressor or activator) and the role of DNA methylation on TF (enhancing or attenuating).

```
results
└── methreg_output
    ├── methylation_matrix_methreg.rds
    ├── expression_matrix_methreg.rds
    ├── triplet_nearest_gene_interactions.tsv
    └── triplet_nearest_gene_interactions_stratified_model.tsv
```

Plots of top 5 triplets per cluster:

```
# n is the cluster number
plots
└── methreg_output
    └── triplet_nearest_gene_interactions_top5_cluster_{n}.pdf
```

### Pathway enrichment of functional CpG sites and differentially expressed target genes

`03-methreg_fgsea_analysis.R`: The purpose of this script is to perform 

1. Pathway enrichment using `missMethyl::gsameth` on the function CpG sites identified by step 2 i.e. `02-methreg_analysis.R`. Only those interactions are used where the the there is a 1) direct effect of DNAm and 2) synergistic effect of DNAm + TF on gene expression.
2. Pathway enrichment using `fgsea` on the target genes of those CpG sites. Only those target genes are used that are found to be differentially expressed using DESeq2 per the `../dge_pathway_analysis` module.

#### Inputs

```
# triplets identified in 02-methreg_analysis.R
results/methreg_output
└── triplet_nearest_gene_interactions.tsv

# differentially expressed genes using DESeq2
../dge_pathway_analysis/results/intNMF/deseq
└── diffexpr_output_per_cluster.tsv

# KEGG Medicus pathway gmt file 
# for differentially methylated CpG sites pathway enrichment
input
└── c2.cp.kegg_medicus.v2023.2.Hs.entrez.gmt

# KEGG Medicus pathway gmt file
# for differentially expressed target genes pathway enrichment
input
└── c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt

# Hallmark pathway gmt file 
# for differentially methylated CpG sites pathway enrichment
input
└── h.all.v2023.2.Hs.entrez.gmt

# Hallmark pathway gmt file
# for differentially expressed target genes pathway enrichment
input
└── h.all.v2023.2.Hs.symbols.gmt
```

#### Outputs

```
results
└── methreg_output
    ├── hallmark
    │   ├── functional_cpg_pathway_enrichment.tsv # enrichment of differentially methylated + functional CpG sites
    │   └── functional_target_gene_pathway_enrichment.tsv # enrichment of differentially expressed target genes 
    └── kegg_medicus
        ├── functional_cpg_pathway_enrichment.tsv
        └── functional_target_gene_pathway_enrichment.tsv
```

### Differential methylated region analysis (DMRcate) + Pathway Enrichment (gsaregion)

`04-dmr_gsaregion_analysis.R`:  The function of this script is to pull the full medulloblastoma methylation m-values dataset and filter the dataset into probes representing promoter region, gene-body (introns and exons) and promoter + gene-body combined. Next, differential region-level methylation analyses is performed with `DMRcate::dmrcate` using a cluster-of-interest vs 'rest' approach, and specifying 'EPIC' array. Finally, `missMethyl::gsaregion` analysis (setting sig.genes = TRUE) is performed on the the resulting differentially methylated regions to determine cluster-specific pathway differences. 

Tests of interest (array.type = 'EPIC'):
a. gene body only
b. promoter only
c. gene body + promoter

#### Inputs

```
# intNMF derived clusters
../intNMF/results/intnmf_clusters.tsv

# methylation annotation file
../../data/v12/infinium.gencode.v39.probe.annotations.tsv.gz

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

### Differential methylated region analysis (DMRcate) + Pathway Enrichment (fgsea)

`05-dmr_fgsea_analysis.R`: The function of this script is to pull the full medulloblastoma methylation m-values dataset and filter the dataset into probes representing promoter region and gene-body (introns and exons). Next, differential region-level methylation analyses is performed with `DMRcate::dmrcate` using a cluster-of-interest vs 'rest' approach, and specifying 'EPIC' array. Finally, `fgsea::fgsea` analysis is performed on the `genes overlapping` the resulting differentially methylated regions to determine cluster-specific pathway differences. 

Tests of interest:
a. gene body only
b. promoter only
c. gene body + promoter

```
# intNMF derived clusters
../intNMF/results/intnmf_clusters.tsv

# methylation annotation file
../../data/v12/infinium.gencode.v39.probe.annotations.tsv.gz

# methylation m-values subsetted to cohort of interest
../data_preparation/data
└── methyl-m-values.rds
```

TSV file of all significant pathways (`FDR < 0.05`) for each cluster:

```
results
└── dmr_fgsea_output
    ├── hallmark
    │   ├── gene_body_fgsea_output_per_cluster.tsv
    │   ├── genebody_promoter_fgsea_output_per_cluster.tsv
    │   └── promoter_fgsea_output_per_cluster.tsv
    └── reactome
        ├── gene_body_fgsea_output_per_cluster.tsv
        ├── genebody_promoter_fgsea_output_per_cluster.tsv
        └── promoter_fgsea_output_per_cluster.tsv
```

Barplots of top 50 pathways (`FDR < 0.05`) enriched in each cluster:
```
plots
└── dmr_fgsea_output
    ├── hallmark
    │   ├── gene_body_fgsea_pathways.pdf
    │   ├── genebody_promoter_fgsea_pathways.pdf
    │   └── promoter_fgsea_pathways.pdf
    └── reactome
        ├── gene_body_fgsea_pathways.pdf
        ├── genebody_promoter_fgsea_pathways.pdf
        └── promoter_fgsea_pathways.pdf
```

### Integration of MethReg output with DIABLO

`06-methreg_diablo_integration.R`: The function of this script is to utilize triplet information from MethReg analysis (i.e. prioritization of differentially methylated CpG sites) and compute an intersection between DIABLO results (with non-zero loadings) and the functional probes that pass the threshold of `RLM_DNAmGroup:TF_fdr` or `RLM_DNAmGroup_fdr` < 0.05. 

#### Inputs

```
# formatted methylation matrix for MethReg analysis
results/methreg_output
└── methylation_matrix_methreg.rds

# formatted expression matrix for MethReg analysis
results/methreg_output
└── expression_matrix_methreg.rds

# output of functionally relevant triplets
results/methreg_output
└── triplet_nearest_gene_interactions_stratified_model.tsv

# output of DIABLO descriptive analysis
../diablo/results/intnmf/descriptive
└── diablo_MB.rds
```

#### Outputs

```
# tsv file of triplets intersecting with DIABLO
results
└── methreg_output
    └── diablo_integration
        └── comp{component_number}_diablo_methreg_intersection.tsv

# interaction model for probes intersecting with DIABLO
plots
└── methreg_output
    └── diablo_integration
        └── comp{component_number}_diablo_methreg_intersection.pdf
```

### Integration of Limma output with DIABLO

`06-limma_diablo_integration.R`: The function of this script is to utilize differentially expressed CpG sites obtained from `01-limma_analysis.R` and compute an intersection with DIABLO results (with non-zero loadings).

#### Inputs

```
# output of differentially expressed CpG sites
results/limma_output
└── genebody_promoter_diffexpr_probes_per_cluster.tsv

# output of DIABLO descriptive analysis
../diablo/results/intnmf/descriptive
└── diablo_MB.rds
```

#### Outputs

```
# tsv file of differential CpGs intersecting with DIABLO
results
└── limma_output
    └── diablo_integration
        └── comp{component_number}_diablo_limma_intersection.tsv
```
