# Workflow Diagram

```mermaid
graph TD
    %% Define main input nodes
    histology[Histology Data]
    counts[RNA-seq Counts]
    methyl[Methylation Beta Values]
    splice[Splicing PSI Values]
    gtf[Gene Annotations]
    
    %% Define parameters
    params[Parameters<br>num_features<br>max_k]
    
    %% Define main workflow steps
    prep[Data Preparation]
    cluster[IntNMF Clustering]
    dge[Differential Expression]
    methyl_analysis[Methylation Analysis]
    post[Post-Clustering Analysis]
    
    %% Define outputs
    clusters[Cluster Assignments]
    feature_scores[Feature Scores]
    plots1[Clustering Plots]
    
    dge_results[DEG Results]
    pathway_rna[RNA Pathways]
    plots2[Expression Plots]
    
    dms_results[DM Results]
    pathway_methyl[Methylation Pathways]
    plots3[Methylation Plots]
    
    associations[Subtype Associations]
    heatmaps[Heatmaps]
    survival[Survival Plots]
    
    %% Connect inputs to steps
    histology --> prep
    counts --> prep
    methyl --> prep
    splice --> prep
    gtf --> prep
    params --> prep
    
    %% Connect steps in main flow
    prep --> |RNA Data| cluster
    prep --> |Methylation Data| cluster
    prep --> |Splicing Data| cluster
    params --> cluster
    
    %% Branching to parallel analyses
    cluster --> clusters
    cluster --> feature_scores
    cluster --> plots1
    
    clusters --> dge
    gtf --> dge
    counts --> dge
    
    clusters --> methyl_analysis
    methyl --> methyl_analysis
    
    clusters --> post
    histology --> post
    feature_scores --> post
    prep -->|All Data| post
    
    %% Connect to outputs
    dge --> dge_results
    dge --> pathway_rna
    dge --> plots2
    
    methyl_analysis --> dms_results
    methyl_analysis --> pathway_methyl
    methyl_analysis --> plots3
    
    post --> associations
    post --> heatmaps
    post --> survival
    
    %% Styling
    classDef input fill:#d1c4e9,stroke:#7e57c2,stroke-width:2px
    classDef process fill:#bbdefb,stroke:#1976d2,stroke-width:2px
    classDef output fill:#c8e6c9,stroke:#388e3c,stroke-width:2px
    classDef param fill:#ffecb3,stroke:#ffa000,stroke-width:2px
    
    class histology,counts,methyl,splice,gtf input
    class prep,cluster,dge,methyl_analysis,post process
    class clusters,feature_scores,plots1,dge_results,pathway_rna,plots2,dms_results,pathway_methyl,plots3,associations,heatmaps,survival output
    class params param
```

## Step-by-Step Workflow Description

### 1. Data Preparation
- **Inputs**: Histology data, RNA-seq counts, methylation beta values, splicing PSI values, gene annotations
- **Parameters**: Number of features to select
- **Process**: Filters data to select samples of interest, identifies most variable features, and performs normalization
- **Outputs**: Filtered and transformed data matrices

### 2. IntNMF Clustering
- **Inputs**: Preprocessed RNA, methylation, and splicing data
- **Parameters**: Maximum number of clusters (k)
- **Process**: Performs integrative non-negative matrix factorization to identify multi-omic clusters
- **Outputs**: Cluster assignments, feature scores for each modality, consensus and silhouette plots

### 3. Differential Gene Expression
- **Inputs**: Cluster assignments, RNA-seq counts, gene annotations
- **Process**: Identifies differentially expressed genes between clusters using DESeq2
- **Outputs**: DEG results, pathway enrichment results, gene expression plots

### 4. Methylation Analysis
- **Inputs**: Cluster assignments, methylation data
- **Process**: Identifies differentially methylated regions and performs methylation-specific pathway analysis
- **Outputs**: DM results, methylation pathway enrichment, methylation plots

### 5. Post-Clustering Analysis
- **Inputs**: Cluster assignments, feature scores, all data modalities, histology data
- **Process**: Performs statistical comparisons with known subtypes, creates visualizations, generates survival plots
- **Outputs**: Association statistics, heatmaps, bubble plots, survival curves, Sankey diagrams

## Data Flow Patterns

The workflow demonstrates several key data flow patterns:

1. **Linear Pipeline**: Data preparation → clustering
2. **Fan-out**: Clustering results feed into three parallel analysis streams
3. **Integration**: Multiple data types are integrated during clustering
4. **Reference Data**: Gene annotations are used at multiple steps