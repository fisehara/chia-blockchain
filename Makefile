DOCKER_REGISTRY=harafise
APP_NAME=chia_blockchain_containered_test
VERSION=alpha-1.3 
BUILD_BRANCH=alpha-1.3
UPNP_PORT=8444
UPNP_PORT_1=8555
UI_PORT=8222

build-docker: ## Build the docker container from dist
	docker build --rm --build-arg BUILD_BRANCH=$(BUILD_BRANCH) -f Dockerfile -t $(DOCKER_REGISTRY)/$(APP_NAME):$(VERSION) .

deploy-docker: ## Build the docker container from dist
	docker push $(DOCKER_REGISTRY)/$(APP_NAME):$(VERSION)

generate-keys:
	docker run --network host $(DOCKER_REGISTRY)/$(APP_NAME):$(VERSION) /bin/bash -c "python -m scripts.regenerate_keys"

run-farmer: ## Run docker image local
	docker run --network host $(DOCKER_REGISTRY)/$(APP_NAME):$(VERSION) /bin/bash -c ". .venv/bin/activate && sh ./scripts/run_farming.sh"

run-full-node: ## Run docker image local
	docker run --network host $(DOCKER_REGISTRY)/$(APP_NAME):$(VERSION) /bin/bash -c ". .venv/bin/activate && sh ./scripts/run_farming.sh"
