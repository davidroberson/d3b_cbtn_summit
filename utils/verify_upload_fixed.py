#!/usr/bin/env python3

import sevenbridges as sb
import sys

# Configuration
PROJECT_ID = 'david.roberson/copy-of-cbtn-summit-2024'
AUTH_TOKEN = 'fa48ce8571484932b25c3d56511b2f95'
FILENAME = 'multi_modal_clustering_workflow_refactored.cwl'

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

# Look for the file
print(f"Looking for file {FILENAME} in project...")
try:
    files = list(api.files.query(project=project.id, names=[FILENAME]))
    if files:
        print(f"Found file: {files[0].name} (ID: {files[0].id})")
        print(f"File size: {files[0].size} bytes")
        print(f"Created on: {files[0].created_on}")
        print(f"Modified on: {files[0].modified_on}")
        print(f"Successful upload confirmed!")
        print(f"\nFile can be accessed at: https://cavatica.sbgenomics.com/u/{PROJECT_ID}/files/{files[0].id}")
    else:
        print(f"File {FILENAME} not found in the project.")
        print("The upload may have failed or is still in progress.")
except Exception as e:
    print(f"Error searching for file: {str(e)}")
    sys.exit(1)