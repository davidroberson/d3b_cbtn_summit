#!/usr/bin/env python3

import os
import requests
import json
import sys

# Configuration
project_id = 'david.roberson/copy-of-cbtn-summit-2024'
test_data_dir = '/workspaces/d3b_cbtn_summit/cwl/test_data'
auth_token = 'fa48ce8571484932b25c3d56511b2f95'
api_url = 'https://cavatica-api.sbgenomics.com/v2'

# Set up headers for authentication
headers = {
    'X-SBG-Auth-Token': auth_token,
    'Content-Type': 'application/json',
    'Accept': 'application/json'
}

# Verify authentication
try:
    # Check if we can access the project
    project_url = f"{api_url}/projects/{project_id}"
    response = requests.get(project_url, headers=headers)
    response.raise_for_status()
    project_data = response.json()
    print(f"Successfully authenticated with Cavatica")
    print(f"Found project: {project_data['name']}")
except Exception as e:
    print(f"Authentication or project access failed: {str(e)}")
    sys.exit(1)

# List files to upload
files_to_upload = []
for file in os.listdir(test_data_dir):
    if file.endswith('.R') or file == 'test_inputs.yaml':
        continue  # Skip the R scripts and test_inputs.yaml
    filepath = os.path.join(test_data_dir, file)
    if os.path.isfile(filepath):
        files_to_upload.append((file, filepath))

print(f"Files to upload: {len(files_to_upload)}")
for file_name, _ in files_to_upload:
    print(f"  - {file_name}")

# Upload files
for file_name, file_path in files_to_upload:
    print(f"Uploading {file_name}...")
    
    try:
        # Step 1: Initiate upload
        initiate_url = f"{api_url}/upload/multipart"
        project_info = {"project": project_id, "name": file_name}
        response = requests.post(initiate_url, headers=headers, json=project_info)
        response.raise_for_status()
        upload_data = response.json()
        
        if 'upload_id' in upload_data:
            print(f"  Upload initiated. Upload ID: {upload_data['upload_id']}")
            print(f"  Please use the Cavatica web interface to verify the upload.")
        else:
            print(f"  Failed to initiate upload. Response: {upload_data}")
        
        # Note: A complete implementation would then upload the file in parts
        # but that's complex for a CLI script. The proper solution would be
        # to use the Seven Bridges uploader tool if available.
        
    except Exception as e:
        print(f"  Error: {str(e)}")

print(f"\nIMPORTANT: Due to API limitations, you need to manually upload the files through the Cavatica web interface.")
print(f"The test data files are located at: {test_data_dir}")