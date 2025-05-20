#!/usr/bin/env python3

import os
import sevenbridges as sb
from sevenbridges.http.error_handlers import rate_limit_sleeper, maintenance_sleeper
import sys

# Configuration
project_name = 'david.roberson/copy-of-cbtn-summit-2024'
test_data_dir = '/workspaces/d3b_cbtn_summit/cwl/test_data'

# Get authentication token from user
auth_token = input("Enter your Cavatica authentication token: ")

# Create API client with the provided token
try:
    api = sb.Api(url='https://cavatica-api.sbgenomics.com/v2', 
                token=auth_token,
                error_handlers=[rate_limit_sleeper, maintenance_sleeper])
    print(f"Successfully authenticated with Cavatica")
except Exception as e:
    print(f"Authentication failed: {str(e)}")
    sys.exit(1)

# Find the project
try:
    project = api.projects.get(id=project_name)
    print(f"Found project: {project.name} (id: {project.id})")
except Exception as e:
    print(f"Error finding project: {str(e)}")
    sys.exit(1)

# List files to upload
files_to_upload = []
for file in os.listdir(test_data_dir):
    if file.endswith('.R') or file == 'test_inputs.yaml':
        continue  # Skip the R scripts and test_inputs.yaml
    filepath = os.path.join(test_data_dir, file)
    if os.path.isfile(filepath):
        files_to_upload.append(filepath)

print(f"Files to upload: {len(files_to_upload)}")
for file in files_to_upload:
    print(f"  - {os.path.basename(file)}")

# Upload files
successful_uploads = 0
failed_uploads = 0

for file_path in files_to_upload:
    file_name = os.path.basename(file_path)
    print(f"Uploading {file_name}...")
    
    try:
        upload = api.files.upload(
            project=project.id,
            path=file_path,
            file_name=file_name,
            overwrite=True  # Overwrite if file exists
        )
        print(f"  Success! File ID: {upload.id}")
        successful_uploads += 1
    except Exception as e:
        print(f"  Failed: {str(e)}")
        failed_uploads += 1

print(f"\nUpload Summary:")
print(f"  Successful: {successful_uploads}")
print(f"  Failed: {failed_uploads}")
print(f"  Total: {len(files_to_upload)}")