### Author: Komal S. Rathi, Adam Kraya

### Purpose

The purpose of this module is to perform Differential Gene Expression and Pathway Analysis on gene expression data.

### Run analysis

To run the full analysis, use the bash script as follows:

```
bash run_analysis.sh
```
***
`01-dge_analysis_deseq.R`:  The function of this script is to perform differential gene expression using [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) R package using the `cluster vs rest` approach. 

#### Inputs

```
../../data
├── c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt # KEGG MEDICUS gmt file
└── gencode.v39.primary_assembly.annotation.gtf.gz # gencode v39

# cohort specific files
../data_preparation/data
└── gene-counts-rsem-expected_count-collapsed.rds

# intNMF inputs but any other clustering method can be used
../intNMF
└── results/intnmf_clusters.tsv # intNMF clusters
```

#### Outputs

For all the `cluster-vs-rest` comparisons, the DESeq2 output with FDR adjusted p-value < 0.05 is saved under `diffexpr_output_per_cluster.tsv` and the pathway enrichment output using `clusterProfiler::GSEA` with FDR adjusted p-value < 0.05 is saved under `_gsea.tsv`.

```
# DESeq2 output
results/intNMF/deseq
├── diffexpr_output_per_cluster.tsv
├── hallmark
│   └── cluster_{n}_vs_rest_gsea.tsv
├── kegg_medicus
│   └── cluster_{n}_vs_rest_gsea.tsv
└── reactome
    └── cluster_{n}_vs_rest_gsea.tsv
```

For each `cluster-vs-rest` enrichment, a barplot of top 10 upregulated and top 10 downregulated (i.e. maximum of 20 pathways) identified by `clusterProfiler::GSEA` utilizing `KEGG MEDICUS` and `REACTOME` pathways at `FDR < 0.05` is generated under `*_gsea_barplot.pdf`, dotplot is generated under `*_gsea_dotplot.pdf` and network under `*_gsea_cnet.pdf`. 

```
# here n is cluster number of interest
# DESeq2 output
plots/intNMF/deseq
├── hallmark
│   ├── cluster_{n}_vs_rest_gsea_barplot.pdf
│   ├── cluster_{n}_vs_rest_gsea_cnet.pdf
│   └── cluster_{n}_vs_rest_gsea_dotplot.pdf
├── kegg_medicus
│   ├── cluster_{n}_vs_rest_gsea_barplot.pdf
│   ├── cluster_{n}_vs_rest_gsea_cnet.pdf
│   └── cluster_{n}_vs_rest_gsea_dotplot.pdf
└── reactome
    ├── cluster_{n}_vs_rest_gsea_barplot.pdf
    ├── cluster_{n}_vs_rest_gsea_cnet.pdf
    └── cluster_{n}_vs_rest_gsea_dotplot.pdf
```

***
`01-dge_analysis_noiseq.R`: The function of this script is to perform differential gene expression using [NOISeq](https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html) R package using the `cluster vs rest` approach.

#### Inputs

```
../../data
├── c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt # KEGG MEDICUS gmt file
└── gencode.v39.primary_assembly.annotation.gtf.gz # gencode v39

# cohort specific files
../data_preparation/data
└── gene-counts-rsem-expected_count-collapsed.rds

# intNMF inputs but any other clustering method can be used
../intNMF
└── results/intnmf_clusters.tsv # intNMF clusters
```

#### Outputs

For all the `cluster-vs-rest` comparisons, the NOISeq output with FDR adjusted p-value < 0.05 is saved under `diffexpr_output_per_cluster.tsv` and the pathway enrichment output using `clusterProfiler::GSEA` with FDR adjusted p-value < 0.05 is saved under `_gsea.tsv`.


```
# NOISeq output
results/intNMF/noiseq
├── diffexpr_output_per_cluster.tsv
├── hallmark
│   └── cluster_{n}_vs_rest_gsea.tsv
├── kegg_medicus
│   └── cluster_{n}_vs_rest_gsea.tsv
└── reactome
    └── cluster_{n}_vs_rest_gsea.tsv
```

For each `cluster-vs-rest` enrichment, a barplot of top 10 upregulated and top 10 downregulated (i.e. maximum of 20 pathways) identified by `clusterProfiler::GSEA` utilizing `KEGG MEDICUS` and `REACTOME` pathways at `FDR < 0.05` is generated under `*_gsea_barplot.pdf`, dotplot is generated under `*_gsea_dotplot.pdf` and network under `*_gsea_cnet.pdf`. 


When using NOISeq for differential expression analysis, `noiseq_pca.pdf` is generated which has the PCA plot of input matrix before and after NOISeq batch correction. 

```
# here n is cluster number of interest
# NOISeq output
plots/intNMF/noiseq
├── hallmark
│   ├── cluster_{n}_vs_rest_gsea_barplot.pdf
│   ├── cluster_{n}_vs_rest_gsea_cnet.pdf
│   └── cluster_{n}_vs_rest_gsea_dotplot.pdf
├── kegg_medicus
│   ├── cluster_{n}_vs_rest_gsea_barplot.pdf
│   ├── cluster_{n}_vs_rest_gsea_cnet.pdf
│   └── cluster_{n}_vs_rest_gsea_dotplot.pdf
├── noiseq_pca.pdf
└── reactome
    ├── cluster_{n}_vs_rest_gsea_barplot.pdf
    ├── cluster_{n}_vs_rest_gsea_cnet.pdf
    └── cluster_{n}_vs_rest_gsea_dotplot.pdf
```
***

`02-dge_diablo_integration.R`: Purpose of this script is to use DGEs from `DESeq2` or `NOISeq` and compute an intersection with the output of `DIABLO` descriptive analysis.

#### Inputs:

```
# DIABLO descriptive output
../diablo/results/intnmf/descriptive
└── diablo_MB.rds

# DESeq2 output
results/intNMF/deseq
└── diffexpr_output_per_cluster.tsv

# or NOISeq output
results/intNMF/noiseq
└── diffexpr_output_per_cluster.tsv
```

#### Outputs:

```
# DESeq2
results/intNMF/deseq/diablo_integration
└── comp{component_number}_diablo_dge_intersection.tsv

# NOISeq
results/intNMF/noiseq/diablo_integration
└── comp{component_number}_diablo_dge_intersection.tsv
```
