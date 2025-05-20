#!/usr/bin/env python3

import os
import sevenbridges as sb
from sevenbridges.http.error_handlers import rate_limit_sleeper, maintenance_sleeper
import sys
import time

# Configuration
project_name = 'david.roberson/copy-of-cbtn-summit-2024'
test_data_dir = '/workspaces/d3b_cbtn_summit/cwl/test_data'
auth_token = 'fa48ce8571484932b25c3d56511b2f95'  # Using previously entered token

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
        # The upload process returns a File object when completed
        file_size = os.path.getsize(file_path)
        
        # For larger files, we'll upload in chunks
        if file_size > 5 * 1024 * 1024:  # 5MB
            print(f"  Large file detected ({file_size/1024/1024:.2f} MB). Using multipart upload...")
            upload = api.files.upload(
                project=project.id,
                path=file_path,
                file_name=file_name,
                overwrite=True,
                part_size=5242880  # 5MB chunks
            )
        else:
            # For smaller files
            upload = api.files.upload(
                project=project.id,
                path=file_path,
                file_name=file_name,
                overwrite=True
            )
            
        # Upload is actually a generator that we need to consume
        # The last yielded value should be the completed file
        completed_file = None
        for f in upload:
            completed_file = f
            # Wait a bit to ensure the file is processed
            time.sleep(1)
        
        if completed_file:
            print(f"  Success! File ID: {completed_file.id}")
            successful_uploads += 1
        else:
            print(f"  Failed: Upload completed but no file was returned")
            failed_uploads += 1
            
    except Exception as e:
        print(f"  Failed: {str(e)}")
        failed_uploads += 1

print(f"\nUpload Summary:")
print(f"  Successful: {successful_uploads}")
print(f"  Failed: {failed_uploads}")
print(f"  Total: {len(files_to_upload)}")