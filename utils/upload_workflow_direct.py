#!/usr/bin/env python3

import os
import sevenbridges as sb
import sys

# Configuration
PROJECT_ID = 'david.roberson/copy-of-cbtn-summit-2024'
AUTH_TOKEN = 'fa48ce8571484932b25c3d56511b2f95'
WORKFLOW_FILE = '/workspaces/d3b_cbtn_summit/cwl/multi_modal_clustering_workflow_refactored.cwl'

# Create API client
api = sb.Api(url='https://cavatica-api.sbgenomics.com/v2', token=AUTH_TOKEN)
print(f"Successfully authenticated with Cavatica")

# Find the project
project = api.projects.get(id=PROJECT_ID)
print(f"Found project: {project.name}")

# Upload the file
print(f"Uploading workflow file to project files...")

# Use the direct upload method
file_path = WORKFLOW_FILE
file_name = os.path.basename(file_path)

try:
    # Upload file using the Files API - using the simple upload endpoint
    upload = api.files.upload(file_path, project, file_name, overwrite=True)
    print(f"Upload initiated for {file_name}")
    print(f"Please check the project files in Cavatica to verify the upload.")
    print(f"The file should be accessible at: https://cavatica.sbgenomics.com/u/{PROJECT_ID}/files/")
except Exception as e:
    print(f"Upload failed: {str(e)}")