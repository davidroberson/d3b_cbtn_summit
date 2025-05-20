#!/bin/bash

# Configuration
PROJECT_ID="david.roberson/copy-of-cbtn-summit-2024"
TEST_DATA_DIR="/workspaces/d3b_cbtn_summit/cwl/test_data"
AUTH_TOKEN="fa48ce8571484932b25c3d56511b2f95"
API_URL="https://cavatica-api.sbgenomics.com/v2"

# Verify authentication and project access
echo "Verifying authentication and project access..."
PROJECT_INFO=$(curl -s -H "X-SBG-Auth-Token: $AUTH_TOKEN" -H "Accept: application/json" "$API_URL/projects/$PROJECT_ID")

if [[ $PROJECT_INFO == *"name"* ]]; then
    PROJECT_NAME=$(echo $PROJECT_INFO | grep -o '"name":"[^"]*"' | cut -d'"' -f4)
    echo "Successfully authenticated with Cavatica"
    echo "Found project: $PROJECT_NAME"
else
    echo "Failed to access project. Response: $PROJECT_INFO"
    exit 1
fi

# List files to upload
echo "Files to upload from $TEST_DATA_DIR:"
files_to_upload=()
for file in "$TEST_DATA_DIR"/*; do
    filename=$(basename "$file")
    # Skip R scripts and YAML files
    if [[ "$filename" != *.R && "$filename" != test_inputs.yaml && -f "$file" ]]; then
        files_to_upload+=("$file")
        echo "  - $filename"
    fi
done

echo ""
echo "IMPORTANT: Due to API limitations, we recommend manually uploading the files through the Cavatica web interface."
echo "The test data files are located at: $TEST_DATA_DIR"
echo ""
echo "To manually upload these files:"
echo "1. Log in to Cavatica at https://cavatica.sbgenomics.com/"
echo "2. Navigate to the project '$PROJECT_NAME'"
echo "3. Click on the 'Files' tab"
echo "4. Click on the 'Upload' button"
echo "5. Select all the test data files and upload them"

# Attempt to initiate upload for each file
echo ""
echo "Attempting to initiate uploads via API anyway (this may fail)..."

for file_path in "${files_to_upload[@]}"; do
    filename=$(basename "$file_path")
    echo "Initiating upload for $filename..."
    
    # Create upload initiation request
    response=$(curl -s -X POST \
        -H "X-SBG-Auth-Token: $AUTH_TOKEN" \
        -H "Content-Type: application/json" \
        -H "Accept: application/json" \
        -d "{\"project\":\"$PROJECT_ID\",\"name\":\"$filename\"}" \
        "$API_URL/upload/multipart")
    
    # Check if the response contains an upload ID
    if [[ $response == *"upload_id"* ]]; then
        upload_id=$(echo $response | grep -o '"upload_id":"[^"]*"' | cut -d'"' -f4)
        echo "  Upload initiated. Upload ID: $upload_id"
    else
        echo "  Failed to initiate upload. Response: $response"
    fi
done