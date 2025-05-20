#!/bin/bash

# Configuration
PROJECT_ID="david.roberson/copy-of-cbtn-summit-2024"
AUTH_TOKEN="fa48ce8571484932b25c3d56511b2f95"
WORKFLOW_FILE="/workspaces/d3b_cbtn_summit/cwl/multi_modal_clustering_workflow_refactored.cwl"
FILE_NAME=$(basename "$WORKFLOW_FILE")

echo "Uploading workflow file $FILE_NAME to CAVATICA..."

# Step 1: Initiate a multipart upload
response=$(curl -s -X POST \
  -H "X-SBG-Auth-Token: $AUTH_TOKEN" \
  -H "Content-Type: application/json" \
  -H "Accept: application/json" \
  -d "{\"project\":\"$PROJECT_ID\",\"name\":\"$FILE_NAME\"}" \
  "https://cavatica-api.sbgenomics.com/v2/upload/multipart")

# Check for errors
if [[ $response == *"message"* && $response == *"error"* ]]; then
  echo "Error initiating upload: $response"
  exit 1
fi

# Check for upload ID
if [[ $response == *"upload_id"* ]]; then
  upload_id=$(echo $response | grep -o '"upload_id":"[^"]*"' | cut -d'"' -f4)
  
  if [ -n "$upload_id" ]; then
    echo "Upload initiated with ID: $upload_id"
    echo "Please check CAVATICA to see if the file appears in the project files."
    echo "You can visit: https://cavatica.sbgenomics.com/u/$PROJECT_ID/files"
  else
    echo "Failed to extract upload ID from response: $response"
    exit 1
  fi
else
  echo "Unexpected response: $response"
  exit 1
fi