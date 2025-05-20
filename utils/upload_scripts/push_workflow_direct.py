#!/usr/bin/env python3

import os
import sevenbridges as sb
import sys
import json

# Configuration
PROJECT_ID = 'david.roberson/copy-of-cbtn-summit-2024'
AUTH_TOKEN = 'fa48ce8571484932b25c3d56511b2f95'
WORKFLOW_FILE = '/workspaces/d3b_cbtn_summit/cwl/multi_modal_clustering_workflow_refactored.cwl'
WORKFLOW_NAME = 'Multi-Modal Clustering Workflow (Refactored)'

# Create API client
try:
    api = sb.Api(url='https://cavatica-api.sbgenomics.com/v2', token=AUTH_TOKEN)
    print(f"Successfully authenticated with Cavatica")
except Exception as e:
    print(f"Authentication failed: {str(e)}")
    sys.exit(1)

# Find the project
try:
    project = api.projects.get(id=PROJECT_ID)
    print(f"Found project: {project.name}")
except Exception as e:
    print(f"Error finding project: {str(e)}")
    sys.exit(1)

# Read workflow file
try:
    with open(WORKFLOW_FILE, 'r') as f:
        workflow_content = f.read()
    print(f"Successfully read workflow file: {WORKFLOW_FILE}")
except Exception as e:
    print(f"Error reading workflow file: {str(e)}")
    sys.exit(1)

# Check if workflow already exists
existing_workflow = None
try:
    workflows = list(api.apps.query(project=project.id, name=WORKFLOW_NAME))
    if workflows:
        existing_workflow = workflows[0]
        print(f"Found existing workflow: {existing_workflow.name} (ID: {existing_workflow.id})")
except Exception as e:
    print(f"Error checking for existing workflow: {str(e)}")

# Create or update the workflow
try:
    if existing_workflow:
        # Update existing workflow
        print(f"Updating existing workflow...")
        app = api.apps.update(
            id=existing_workflow.id,
            raw={"app": workflow_content}
        )
        print(f"Successfully updated workflow: {app.name} (ID: {app.id})")
    else:
        # Create new workflow
        print(f"Creating new workflow...")
        app = api.apps.create_app(
            project=project.id,
            name=WORKFLOW_NAME,
            raw={"app": workflow_content}
        )
        print(f"Successfully created workflow: {app.name} (ID: {app.id})")
    
    print(f"\nWorkflow can be accessed at: https://cavatica.sbgenomics.com/u/{PROJECT_ID}/apps/{app.id}")
except Exception as e:
    print(f"Error pushing workflow: {str(e)}")
    sys.exit(1)