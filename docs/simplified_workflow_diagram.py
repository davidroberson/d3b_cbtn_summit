#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch

def create_simplified_diagram(output_file):
    """Create a simplified workflow diagram focused on the key pipeline steps."""
    
    # Setup figure and grid
    fig = plt.figure(figsize=(12, 8), dpi=300)
    gs = gridspec.GridSpec(3, 5, height_ratios=[1, 3, 1])
    
    # Define colors
    colors = {
        'input': '#BBDEFB',       # Light blue
        'prep': '#C8E6C9',        # Light green
        'cluster': '#FFCDD2',     # Light red
        'analysis': '#FFE0B2',    # Light orange
        'output': '#E1BEE7',      # Light purple
        'arrow': '#757575',       # Dark gray
        'env': '#ECEFF1'          # Light gray
    }
    
    # Title
    ax_title = fig.add_subplot(gs[0, 1:4])
    ax_title.text(0.5, 0.5, 'CBTN Multi-Omic Clustering Workflow', 
                 fontsize=18, fontweight='bold', ha='center', va='center')
    ax_title.axis('off')
    
    # Main workflow sections
    sections = [
        {'name': 'INPUT DATA', 'color': colors['input'], 'pos': gs[1, 0],
         'content': "• Histology Data\n• RNA-seq Counts\n• Methylation Values\n• Splicing PSI Values\n• Gene Annotations"},
        
        {'name': 'DATA PREPARATION', 'color': colors['prep'], 'pos': gs[1, 1],
         'content': "1. Filter samples by\n   histology type\n2. Select top variable\n   features\n3. Transform and\n   normalize data"},
        
        {'name': 'INTNMF CLUSTERING', 'color': colors['cluster'], 'pos': gs[1, 2],
         'content': "1. Prepare data matrices\n2. Run IntNMF with\n   multiple K values\n3. Select optimal K\n4. Generate feature\n   importance scores"},
        
        {'name': 'DOWNSTREAM ANALYSIS', 'color': colors['analysis'], 'pos': gs[1, 3],
         'content': "• Differential Gene Expr.\n  (DESeq2 analysis)\n• Methylation Analysis\n  (Limma analysis)\n• Post-clustering\n  - Statistical tests\n  - Visualizations\n  - Survival analysis"},
        
        {'name': 'OUTPUTS', 'color': colors['output'], 'pos': gs[1, 4],
         'content': "• Cluster assignments\n• Gene expression\n  signatures\n• Methylation patterns\n• Survival plots\n• Pathway enrichment\n  results"}
    ]
    
    # Draw section boxes
    axes = []
    for section in sections:
        ax = fig.add_subplot(section['pos'])
        ax.set_facecolor(section['color'])
        
        # Add a border effect
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('#424242')
            spine.set_linewidth(1.5)
        
        # Title on top of box
        ax.text(0.5, 0.95, section['name'], fontsize=11, 
                fontweight='bold', ha='center', va='top',
                bbox=dict(facecolor='white', alpha=0.8, pad=3, 
                          edgecolor='#424242', boxstyle='round,pad=0.5'))
        
        # Content in box
        ax.text(0.5, 0.5, section['content'], fontsize=9,
                ha='center', va='center', linespacing=1.5)
        
        ax.set_xticks([])
        ax.set_yticks([])
        axes.append(ax)
    
    # Add workflow arrows between sections
    for i in range(len(axes) - 1):
        fig.add_artist(FancyArrowPatch(
            posA=(axes[i].get_position().x1, axes[i].get_position().y0 + axes[i].get_position().height/2),
            posB=(axes[i+1].get_position().x0, axes[i+1].get_position().y0 + axes[i+1].get_position().height/2),
            arrowstyle='-|>', color=colors['arrow'],
            connectionstyle='arc3,rad=0.0', linewidth=2,
            mutation_scale=20, zorder=1000
        ))
    
    # Environment box at bottom
    ax_env = fig.add_subplot(gs[2, :])
    ax_env.set_facecolor(colors['env'])
    for spine in ax_env.spines.values():
        spine.set_visible(True)
        spine.set_color('#424242')
        
    # Environment content
    env_text = (
        "WORKFLOW EXECUTION ENVIRONMENT\n"
        "Container: bioconductor/bioconductor_docker:RELEASE_3_18 | CWL Execution: cwltool runtime\n"
        "Core: tidyverse, IntNMF, survival | Bioconductor: DESeq2, limma, ComplexHeatmap"
    )
    
    ax_env.text(0.5, 0.5, env_text, fontsize=9,
            ha='center', va='center', linespacing=1.5,
            bbox=dict(facecolor='white', alpha=0.8, pad=6, 
                      edgecolor='#424242', boxstyle='round,pad=0.5'))
    
    ax_env.set_xticks([])
    ax_env.set_yticks([])
    
    # Description at the very bottom
    fig.text(0.5, 0.02, 
             "A workflow for multi-modal clustering analysis of cancer genomics data using Common Workflow Language (CWL)",
             ha='center', fontsize=9, style='italic')
    
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.1, hspace=0.3)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Simplified workflow diagram saved to {output_file}")

if __name__ == "__main__":
    import os
    
    # Setup file paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = os.path.join(script_dir, "simplified_workflow_diagram.png")
    
    # Create the diagram
    create_simplified_diagram(output_file)