
### Author: Komal S. Rathi

### Purpose

The purpose of this module is:
1. Subset input data matrices (RNA, CNV, SNV, Methylation and Splicing data) to short histology of interest.
2. Format PSI matrix to a format accepted by PEGASAS.
3. Filter/Transform data matrices which can be used as inputs for multi-modal clustering packages. 

### Data version

- OPC v15 for histologies (molecular subtype etc), RNA, Methylation, Splicing

### Run Analysis
```
# run full analysis
bash run_analysis.sh
```

### Description of scripts
***

`01-multi-modal-clustering-prepare-data.R` : The purpose of this script is to prepare input files from all available modalities and use them as input for the IntNMF clustering analysis.

#### Feature selection

Features were selected from OpenPedCan-analysis v15 datasets using the following filters:

1) **RNA**: 

- The expected counts dataset was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable protein coding genes` followed by `Rank transformation`.

2) **Methylation**:

- Methylation beta-values matrix was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable probes`.

3) **Splicing**:

- Splice matrix was first filtered to Medulloblatoma samples. 
- Features were reduced to `Top 1000 most variable splice variants`.

#### Output

The script resulted in `152 Medulloblastoma samples` that have all 3 data modalities available. The filtered/transformed data matrices are written out to individual .tsv files. Additionally a mapping between `Kids_First_Biospecimen_ID` identifiers from each modality and `sample_id` is written out to  `samples_map.tsv`.

```
results
├── cnv_data.tsv # cnv data 
├── methyl_data.tsv # methylation data 
├── norm_counts.tsv # expression data 
├── snv_data.tsv # snv data 
├── splice_data.tsv # splice data 
└── samples_map.tsv # biospecimens + cohort identifiers for samples used for each modality 
```
