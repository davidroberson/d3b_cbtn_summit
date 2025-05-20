#!/usr/bin/env python3

import os
import sevenbridges as sb
import sys

# Configuration
PROJECT_ID = 'david.roberson/copy-of-cbtn-summit-2024'
AUTH_TOKEN = 'fa48ce8571484932b25c3d56511b2f95'  # Using previously entered token
WORKFLOW_FILE = '/workspaces/d3b_cbtn_summit/cwl/multi_modal_clustering_workflow_refactored.cwl'

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

# Check if file exists
if not os.path.isfile(WORKFLOW_FILE):
    print(f"Error: Workflow file {WORKFLOW_FILE} does not exist")
    sys.exit(1)

# Upload the file
filename = os.path.basename(WORKFLOW_FILE)
print(f"Uploading {filename} to CAVATICA...")

try:
    # Start the upload
    api.files.upload(
        path=WORKFLOW_FILE,
        project=project.id,
        file_name=filename,
        overwrite=True
    )
    print(f"Upload initiated. Please check CAVATICA to verify the file was uploaded successfully.")
    print(f"Note: The file will be overwritten if it already exists.")
    
except Exception as e:
    print(f"Upload failed: {str(e)}")
    sys.exit(1)

print(f"\nTo check the uploaded file, visit:")
print(f"https://cavatica.sbgenomics.com/u/{PROJECT_ID}/files")