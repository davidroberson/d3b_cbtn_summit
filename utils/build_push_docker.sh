#!/bin/bash

# Ensure the script stops on first error
set -e

# Define variables
DOCKER_REPO="pgc-images.sbgenomics.com/david.roberson/cbtn-multiomic-clustering"
VERSION="v1.0.0"
TAG="${DOCKER_REPO}:${VERSION}"

echo "Building Docker image: ${TAG}"

# Build the Docker image
docker build -t ${TAG} .

echo "Docker image built successfully"

# Optional: List built images to confirm
docker images | grep cbtn-multiomic-clustering

# Login to PGC Docker registry (requires your CAVATICA credentials)
echo "Please login to the PGC Docker registry with your CAVATICA credentials"
echo "Use: docker login pgc-images.sbgenomics.com -u <USERNAME> -p <YOUR-AUTH-TOKEN>"

# Push the image to the repository
echo "Pushing Docker image to ${DOCKER_REPO}"
docker push ${TAG}

echo "Docker image pushed successfully"