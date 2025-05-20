#!/usr/bin/env python3

import os
import sevenbridges as sb
import time
import sys

# Configuration
project_name = 'david.roberson/copy-of-cbtn-summit-2024'
test_data_dir = '/workspaces/d3b_cbtn_summit/cwl/test_data'
auth_token = 'fa48ce8571484932b25c3d56511b2f95'  # Using previously entered token

# Create API client with the provided token
try:
    api = sb.Api(url='https://cavatica-api.sbgenomics.com/v2', token=auth_token)
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
        # Start the upload
        upload_job = api.files.upload(project.id, file_path)
        
        # Wait for upload to complete
        print(f"  Upload initiated. Waiting for completion...")
        
        # Check if the file already exists in the project
        existing_files = api.files.query(project=project.id, names=[file_name])
        
        if len(list(existing_files)) > 0:
            print(f"  Success! File already exists in the project.")
            successful_uploads += 1
            continue
        
        # Wait a bit and check if the file appears
        max_attempts = 10
        for attempt in range(max_attempts):
            time.sleep(5)  # wait 5 seconds between checks
            files = list(api.files.query(project=project.id, names=[file_name]))
            if files:
                print(f"  Success! File ID: {files[0].id}")
                successful_uploads += 1
                break
            else:
                print(f"  Waiting for upload to complete... Attempt {attempt+1}/{max_attempts}")
                
        if attempt == max_attempts - 1 and not files:
            print(f"  Failed: Upload timed out")
            failed_uploads += 1
            
    except Exception as e:
        print(f"  Failed: {str(e)}")
        failed_uploads += 1

print(f"\nUpload Summary:")
print(f"  Successful: {successful_uploads}")
print(f"  Failed: {failed_uploads}")
print(f"  Total: {len(files_to_upload)}")