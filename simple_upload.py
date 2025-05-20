#!/usr/bin/env python3

import os
import sevenbridges as sb
import sys

# Configuration
project_id = 'david.roberson/copy-of-cbtn-summit-2024'
test_data_dir = '/workspaces/d3b_cbtn_summit/cwl/test_data'
auth_token = 'fa48ce8571484932b25c3d56511b2f95'

# Create API client
api = sb.Api(url='https://cavatica-api.sbgenomics.com/v2', token=auth_token)
print(f"Successfully authenticated with Cavatica")

# Get the project
project = api.projects.get(id=project_id)
print(f"Found project: {project.name}")

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
        api.files.upload(file_path, project, file_name)
        print(f"  Upload initiated for {file_name}")
    except Exception as e:
        print(f"  Error: {str(e)}")

print(f"\nUploads initiated. Please check the project on Cavatica to verify the files were uploaded successfully.")