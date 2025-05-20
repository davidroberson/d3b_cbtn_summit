#!/usr/bin/env python3

import os
import json
import subprocess

# Create credentials directory if it doesn't exist
home_dir = os.path.expanduser("~")
credentials_dir = os.path.join(home_dir, ".sevenbridges")
os.makedirs(credentials_dir, exist_ok=True)

# Path to the credentials file
credentials_file = os.path.join(credentials_dir, "credentials")

# Get authentication token from user
token = input("Enter your Cavatica authentication token: ")

# Create the credentials configuration
credentials_config = {
    "default": {
        "api_endpoint": "https://cavatica-api.sbgenomics.com/v2",
        "auth_token": token
    }
}

# Write credentials to file
with open(credentials_file, 'w') as f:
    json.dump(credentials_config, f, indent=4)

# Set proper permissions
subprocess.run(['chmod', '600', credentials_file])

print(f"Credentials file created at: {credentials_file}")