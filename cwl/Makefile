# Makefile for CBTN Multi-Omic Clustering Workflow

.PHONY: validate test docker push clean

# Default CWL runner
CWL_RUNNER = cwltool
DOCKER_REGISTRY = pgc-images.sbgenomics.com/david.roberson
DOCKER_IMAGE = cbtn-multiomic-clustering
DOCKER_TAG = v1.0.0

# Paths
WORKFLOW = multi_modal_clustering_workflow.cwl
TEST_INPUTS = test_data/test_inputs.yaml
TOOLS_DIR = tools
CONTAINER_DIR = containers

# Validate all CWL files
validate:
	@echo "Validating workflow..."
	$(CWL_RUNNER) --validate $(WORKFLOW)
	@echo "Validating tools..."
	for tool in $(TOOLS_DIR)/*.cwl; do \
		echo "Validating $$tool"; \
		$(CWL_RUNNER) --validate $$tool; \
	done

# Generate test data
generate-test-data:
	@echo "Generating test data..."
	docker run --rm -v $(shell pwd):/data -w /data \
		bioconductor/bioconductor_docker:RELEASE_3_18 \
		Rscript test_data/create_test_data_simple.R

# Run workflow with test data
test: validate
	@echo "Running workflow with test data..."
	$(CWL_RUNNER) $(WORKFLOW) $(TEST_INPUTS)

# Run individual workflow steps for testing
test-data-preparation:
	$(CWL_RUNNER) $(TOOLS_DIR)/data_preparation.cwl $(TEST_INPUTS)

test-clustering:
	$(CWL_RUNNER) $(TOOLS_DIR)/integrative_nmf.cwl clustering_inputs.yaml

# Build Docker container
docker:
	@echo "Building Docker image..."
	cd $(CONTAINER_DIR) && \
	docker build -f Dockerfile.simple -t $(DOCKER_REGISTRY)/$(DOCKER_IMAGE):$(DOCKER_TAG) .

# Push Docker container to registry
push: docker
	@echo "Pushing Docker image to registry..."
	docker push $(DOCKER_REGISTRY)/$(DOCKER_IMAGE):$(DOCKER_TAG)

# Upload workflow to CAVATICA
upload-cavatica:
	@echo "Uploading workflow to CAVATICA..."
	sbpack cavatica david.roberson/copy-of-cbtn-summit-2024/multi-modal-clustering-$(DOCKER_TAG) $(WORKFLOW)

# Clean up temporary files
clean:
	@echo "Cleaning up temporary files..."
	rm -rf tmp/
	find . -name "*.tmp" -delete
	find . -name "*.log" -delete

# Show help information
help:
	@echo "CBTN Multi-Omic Clustering Workflow"
	@echo ""
	@echo "Available commands:"
	@echo "  make validate              - Validate all CWL files"
	@echo "  make generate-test-data    - Generate test data"
	@echo "  make test                  - Run workflow with test data"
	@echo "  make test-data-preparation - Test data preparation step only"
	@echo "  make test-clustering       - Test clustering step only"
	@echo "  make docker                - Build Docker image"
	@echo "  make push                  - Push Docker image to registry"
	@echo "  make upload-cavatica       - Upload workflow to CAVATICA"
	@echo "  make clean                 - Clean up temporary files"