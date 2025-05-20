# Multi-Omic Clustering Workflow Visualization

This directory contains visualizations of the CBTN Multi-Omic Clustering Workflow.

## Workflow Diagrams

Two workflow diagrams are available:

1. **Simplified Workflow Diagram** - A clear, high-level representation of the pipeline showing the main processing steps and data flow.
   ![Simplified Workflow Diagram](simplified_workflow_diagram.png)

2. **Detailed Workflow Diagram** - A comprehensive view based directly on the CWL workflow definition, showing all inputs, steps, and outputs with their interconnections.
   ![Detailed Workflow Diagram](workflow_diagram.png)

3. **Mermaid Diagram** - An interactive diagram using the Mermaid.js format (compatible with GitHub markdown).

```mermaid
graph TD
    subgraph "Main Inputs"
        A1[histology_file] --> A[Data Preparation]
        A2[short_histology: HGAT] --> A
        A3[count_file: RNA-seq] --> A
        A4[methyl_file: Methylation] --> A
        A5[splice_file: Splicing] --> A
        A6[gtf_file: Gene Annotation] --> A
        A7[num_features: 1000] --> A
    end

    subgraph "Data Preparation"
        A --> B1[prepared_rna_data]
        A --> B2[prepared_methyl_data]
        A --> B3[prepared_splice_data]
        A --> B4[sample_mapping]
    end

    B1 --> B[Integrative NMF Clustering]
    B2 --> B
    B3 --> B
    B4 --> B
    A8[max_k: 10] --> B

    subgraph "Clustering Results"
        B --> C1[cluster_assignments]
        B --> C2[feature_scores_rna]
        B --> C3[feature_scores_methyl]
        B --> C4[feature_scores_splice]
        B --> C5[consensus_plot]
        B --> C6[silhouette_plot]
    end

    C1 --> C[Differential Expression]
    A3 --> C
    A6 --> C
    A1 --> C

    C1 --> D[Methylation Analysis]
    A4 --> D
    A1 --> D

    C1 --> E[Post-Clustering Analysis]
    B1 --> E
    B2 --> E
    B3 --> E
    C2 --> E
    C3 --> E
    C4 --> E
    A1 --> E

    subgraph "Expression Results"
        C --> F1[deseq_results]
        C --> F2[expression_pathway_results]
        C --> F3[expression_pathway_plots]
    end

    subgraph "Methylation Results"
        D --> G1[methylation_results]
        D --> G2[methylation_pathway_results]
        D --> G3[methylation_pathway_plots]
    end

    subgraph "Post-Clustering Results"
        E --> H1[comparison_results]
        E --> H2[heatmaps]
        E --> H3[bubble_plots]
        E --> H4[survival_plots]
        E --> H5[sankey_plots]
    end
```

## How to Generate Diagrams

The diagrams are generated using Python scripts in this directory:

```bash
# Generate the detailed diagram
python3 generate_workflow_diagram.py

# Generate the simplified diagram
python3 simplified_workflow_diagram.py
```

### Requirements

To generate the diagrams, you need Python with the following libraries:
- networkx
- matplotlib
- pyyaml

Install them with:
```bash
pip install networkx matplotlib pyyaml
```

## Text-Based Representation

A text-based ASCII representation is also available in [workflow_diagram.txt](workflow_diagram.txt), providing a quick overview of the workflow structure that can be viewed directly in the terminal or text editor.

## Workflow Components

### 1. Data Preparation
- **Purpose**: Filter and transform data for clustering
- **Inputs**: Histology data, gene counts, methylation values, splice data, GTF file
- **Outputs**: Filtered data matrices for RNA, methylation, and splicing
- **Key operations**: 
  - Subset to specific histology (HGAT)
  - Select top N variable features
  - Transform data

### 2. Integrative NMF Clustering
- **Purpose**: Perform multi-modal clustering
- **Inputs**: Processed RNA, methylation, and splicing data matrices
- **Outputs**: Cluster assignments, feature importance scores, quality plots
- **Key operations**:
  - Run IntNMF with different k values
  - Select optimal number of clusters
  - Generate cluster assignments and silhouette scores

### 3. Differential Expression Analysis
- **Purpose**: Identify cluster-specific gene expression patterns
- **Inputs**: RNA data, cluster assignments
- **Outputs**: Differentially expressed genes, pathway enrichment results
- **Key operations**:
  - DESeq2 differential expression
  - GSEA pathway analysis
  - Visualization of enriched pathways

### 4. Methylation Analysis
- **Purpose**: Identify cluster-specific methylation patterns
- **Inputs**: Methylation data, cluster assignments
- **Outputs**: Differentially methylated probes, pathway enrichment results
- **Key operations**:
  - Limma differential methylation
  - Methylation-specific pathway analysis
  - Visualization of enriched pathways

### 5. Post-Clustering Analysis
- **Purpose**: Generate visualizations and statistics for clusters
- **Inputs**: Cluster assignments, feature scores, original data
- **Outputs**: Heatmaps, bubble plots, survival curves, etc.
- **Key operations**:
  - Statistical comparisons between clusters and known subtypes
  - Survival analysis
  - Visualization of cluster features

## Data Flow

1. **Input Data** → **Data Preparation**
   - Raw data matrices are filtered and transformed

2. **Prepared Data** → **Clustering**
   - Filtered matrices are used for integrative clustering

3. **Cluster Results** → **Downstream Analyses**
   - Cluster assignments drive differential expression, methylation, and post-clustering analyses

4. **Final Results**
   - Multiple visualization outputs and statistical results
   - Complete characterization of molecular subtypes

## Container Requirements

The workflow uses a Docker container with the following key R packages:
- IntNMF (clustering)
- DESeq2 (differential expression)
- limma (methylation analysis)
- clusterProfiler (pathway analysis)
- survival/survminer (survival analysis)
- ggplot2/ggpubr (visualization)
- ComplexHeatmap (heatmap generation)