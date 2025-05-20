#!/usr/bin/env python3

import os
import yaml
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import matplotlib.patches as mpatches

def parse_cwl(file_path):
    """Parse a CWL file and extract workflow information."""
    with open(file_path, 'r') as f:
        # Skip the shebang line if present
        first_line = f.readline()
        if not first_line.startswith('cwlVersion'):
            cwl_content = first_line + f.read()
        else:
            cwl_content = first_line + f.read()
            
    # Parse YAML content
    try:
        workflow = yaml.safe_load(cwl_content)
        return workflow
    except yaml.YAMLError as e:
        print(f"Error parsing CWL file: {e}")
        return None

def create_workflow_graph(workflow):
    """Create a directed graph representing the workflow."""
    G = nx.DiGraph()
    
    # Add input nodes
    inputs = {}
    for input_id, input_def in workflow.get('inputs', {}).items():
        node_id = f"input:{input_id}"
        label = input_id
        doc = input_def.get('doc', '')
        if doc:
            label = f"{label}\n({doc.split('.')[0]})"
        inputs[input_id] = node_id
        G.add_node(node_id, type='input', label=label, doc=doc)
    
    # Add step nodes and connect from inputs
    steps = {}
    for step_id, step_def in workflow.get('steps', {}).items():
        node_id = f"step:{step_id}"
        steps[step_id] = node_id
        
        # Get the tool/workflow that this step runs
        run_file = step_def.get('run', '')
        if isinstance(run_file, str) and '/' in run_file:
            run_file = os.path.basename(run_file).replace('.cwl', '')
        else:
            run_file = str(run_file)
            
        label = f"{step_id}\n({run_file})"
        G.add_node(node_id, type='step', label=label, run=run_file)
        
        # Connect inputs to this step
        for input_id, input_source in step_def.get('in', {}).items():
            if isinstance(input_source, dict) and 'source' in input_source:
                input_source = input_source['source']
            
            if isinstance(input_source, str):
                source_parts = input_source.split('/')
                if len(source_parts) == 1:
                    # Direct connection from workflow input
                    if input_source in inputs:
                        G.add_edge(inputs[input_source], node_id, 
                                  label=input_id, type='input-to-step')
                elif len(source_parts) == 2:
                    # Connection from another step
                    source_step, source_output = source_parts
                    if source_step in steps:
                        G.add_edge(steps[source_step], node_id, 
                                  label=f"{source_output} → {input_id}", 
                                  type='step-to-step')
    
    # Add output nodes and connect from steps
    outputs = {}
    for output_id, output_def in workflow.get('outputs', {}).items():
        node_id = f"output:{output_id}"
        outputs[output_id] = node_id
        label = output_id
        G.add_node(node_id, type='output', label=label)
        
        output_source = output_def.get('outputSource', '')
        if output_source:
            source_parts = output_source.split('/')
            if len(source_parts) == 2:
                source_step, source_output = source_parts
                if source_step in steps:
                    G.add_edge(steps[source_step], node_id, 
                              label=source_output, type='step-to-output')
    
    return G

def visualize_workflow(G, output_file):
    """Visualize the workflow graph and save to a file."""
    plt.figure(figsize=(16, 12), dpi=100)
    
    # Define node colors
    node_colors = {
        'input': '#BBDEFB',  # Light blue
        'step': '#C8E6C9',   # Light green
        'output': '#FFCC80'  # Light orange
    }
    
    # Position nodes using hierarchical layout
    pos = nx.multipartite_layout(G, subset_key='type', align='horizontal')
    
    # Adjust vertical positions to avoid overlaps
    types = nx.get_node_attributes(G, 'type')
    type_groups = {'input': [], 'step': [], 'output': []}
    for node, ntype in types.items():
        type_groups[ntype].append(node)
    
    # Calculate vertical spacing for each group
    for ntype, nodes in type_groups.items():
        count = len(nodes)
        if count > 1:
            step = 1.0 / (count - 1) if count > 1 else 0
            for i, node in enumerate(sorted(nodes)):
                pos[node] = (pos[node][0], -0.5 + i * step)
    
    # Draw edges with labels
    edge_labels = nx.get_edge_attributes(G, 'label')
    nx.draw_networkx_edges(G, pos, alpha=0.5, node_size=3000, min_source_margin=20, min_target_margin=20)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)
    
    # Draw nodes with custom colors
    for node_type, color in node_colors.items():
        node_list = [n for n in G.nodes() if G.nodes[n]['type'] == node_type]
        nx.draw_networkx_nodes(G, pos, nodelist=node_list, node_color=color, node_size=3000, alpha=0.8)
    
    # Draw node labels
    node_labels = {n: G.nodes[n]['label'] for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=9, font_weight='bold')
    
    # Create legend
    legend_patches = [
        mpatches.Patch(color=node_colors['input'], label='Input Data'),
        mpatches.Patch(color=node_colors['step'], label='Processing Steps'),
        mpatches.Patch(color=node_colors['output'], label='Output Results')
    ]
    plt.legend(handles=legend_patches, loc='upper center', bbox_to_anchor=(0.5, 0.1))
    
    # Add a title
    plt.title('CBTN Multi-Omic Clustering Workflow Diagram', fontsize=16, pad=20)
    
    # Add description text
    description = (
        "A workflow for multi-modal clustering analysis of cancer genomics data. "
        "Takes RNA-seq expression, methylation, and alternative splicing data "
        "as input to generate integrated clusters and perform downstream analyses."
    )
    plt.figtext(0.5, 0.02, description, ha='center', fontsize=10, 
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.7))
    
    # Remove axis
    plt.axis('off')
    plt.tight_layout()
    
    # Save to file
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Workflow diagram saved to {output_file}")

if __name__ == "__main__":
    # Setup file paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    cwl_dir = os.path.dirname(script_dir)
    cwl_file = os.path.join(cwl_dir, "multi_modal_clustering_workflow.cwl")
    output_file = os.path.join(script_dir, "workflow_diagram.png")
    
    # Parse workflow and create graph
    workflow = parse_cwl(cwl_file)
    if workflow:
        G = create_workflow_graph(workflow)
        visualize_workflow(G, output_file)
    else:
        print("Failed to parse workflow file.")