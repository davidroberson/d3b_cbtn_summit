
### Author: Komal S. Rathi

### Purpose

The purpose of this module is:
1. Subset input data matrices (RNA, CNV, SNV, Methylation and Splicing data) to short histology of interest.
2. Format PSI matrix to a format accepted by PEGASAS.
3. Filter/Transform data matrices which can be used as inputs for multi-modal clustering packages. 

### Data version

- OPC v12 for histologies (molecular subtype etc), RNA, Methylation, Splicing and SNV
- [Release 20230309](https://cavatica.sbgenomics.com/u/d3b-bixu-ops/monthly-release-data/files/#q?path=20230309_release) for the CNV gainloss file (as this is not available in OPC).

### Run Analysis
```
# run full analysis
bash run_analysis.sh
```

### Description of scripts
***

`01-subset-data.R`: Purpose of this script is to subset the input data matrices to a specific short histology in order to create subsetted data matrices. These matrices will be used as input files for all downstream analyses modules. 

#### Output

```
data
├── All.gainloss.txt.gz
├── gene-counts-rsem-expected_count-collapsed.rds
├── gene-expression-rsem-tpm-collapsed.rds
├── histologies.tsv
├── methyl-beta-values.rds
├── methyl-m-values.rds
├── snv-consensus-plus-hotspots.maf.tsv.gz
└── splice_events_pan_cancer_functional_filter.rds
```

***
`02-format-splice-data-pegasas.R`: Purpose of this script is to read RMATs file, , and filter to input functional sites per @naqvia and format it per PEGASAS specifications.

#### Output

```
data
└── splice-events-rmats-functional-sites.tsv.gz
```

***
`03-multi-modal-clustering-prepare-data.R` : The purpose of this script is to prepare input files from all available modalities and use them as input for the IntNMF clustering analysis.

#### Feature selection

Features were selected from OpenPedCan-analysis v12 datasets using the following filters:

1) **CNV**:

- The `gainloss.txt` file was first filtered to Medulloblatoma samples.
- Adjusted copy number was calculated. 
- Neutral sites were removed.
- The dataset was reduced to features corresponding to "Gain", "Amplification" in Oncogenes and "Loss", "Complete Loss" in TSGs (using the comprehensive cancer gene list from `annoFuse`)
- Finally, features with `standard deviation >= 0.9` were chosen for clustering. (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0176278)

2) **SNV**:

- The Consensus MAF dataset was first filtered to Medulloblatoma samples.
- Features were reduced to non-synonymous variant classifications e.g. `Missense_Mutation, Frame_Shift_Del, In_Frame_Ins, Frame_Shift_Ins, Splice_Site, Nonsense_Mutation, In_Frame_Del, Nonstop_Mutation, Translation_Start_Site`. 
- Finally, features with `mutation rate > 0.02` (i.e. filter out features with low mutation rate) were chosen for clustering. (Reference: https://www.bioconductor.org/packages/release/bioc/vignettes/iClusterPlus/inst/doc/iManual.pdf)

3) **RNA**: 

- The expected counts dataset was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable protein coding genes` followed by `Rank transformation`.

4) **Methylation**:

- Methylation beta-values matrix was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable probes`.

5) **Splicing**:

- Splice matrix was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable splice variants`.

#### Output

The script resulted in `152 Medulloblastoma samples` that have all 5 data modalities available. The filtered/transformed data matrices are written out to individual .tsv files. Additionally a mapping between `Kids_First_Biospecimen_ID` identifiers from each modality and `sample_id` is written out to  `samples_map.tsv`.

```
results
├── cnv_data.tsv # cnv data 
├── methyl_data.tsv # methylation data 
├── norm_counts.tsv # expression data 
├── snv_data.tsv # snv data 
├── splice_data.tsv # splice data 
└── samples_map.tsv # biospecimens + cohort identifiers for samples used for each modality 
```
