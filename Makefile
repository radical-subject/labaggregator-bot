THIS_FILE := $(lastword $(MAKEFILE_LIST))
MONGOD_STARTED := $(shell systemctl is-active mongod)
DOCKER_COMPOSE_CMD := docker compose -f docker-compose.yml
DOCKER_COMPOSE_DEV_CMD := docker compose -f docker-compose.yml -f docker-compose.development.yml
DOCKER_COMPOSE_TEST_CMD := docker compose -f docker-compose.yml -f docker-compose.pytest.yml

#  COLORS
GREEN  := $(shell tput -Txterm setaf 2)
YELLOW := $(shell tput -Txterm setaf 3)
WHITE  := $(shell tput -Txterm setaf 7)
RESET  := $(shell tput -Txterm sgr0)

HELP_TARGET_MAX_CHAR_NUM=20

.PHONY: \
	pull_n_up \
	pull \
	mongod_stop \
	help \
	dev \
	env_prod \
	env_test \
	goto_app_src \
	build \
	up \
	start \
	down \
	destroy \
	stop \
	restart \
	logs \
	logs_api \
	ps \
	login-timescale \
	login-api \
	db-shell \
	prod_remote \
	mongodb_shell \
	annihilate \
	bot_shell

#  Stops mongod service if it is running
mongod_stop:
	@if [ "$(MONGOD_STARTED)" == "active" ]; then \
		echo "mongod is running; stopping..."; \
		sudo service mongod stop; \
	else \
		echo "mongod not running"; \
	fi


#  Show help
help:
	@echo ''
	@echo 'Usage:'
	@echo '  ${YELLOW}make${RESET} ${GREEN}<target>${RESET}'
	@echo ''
	@echo 'Targets:'
	@awk '/^[a-zA-Z\-\_0-9]+:/ { \
		helpMessage = match(lastLine, /^# (.*)/); \
		if (helpMessage) { \
			helpCommand = substr($$1, 0, index($$1, ":")-1); \
			helpMessage = substr(lastLine, RSTART + 3, RLENGTH); \
			printf "  ${YELLOW}%-$(HELP_TARGET_MAX_CHAR_NUM)s${RESET} ${GREEN}%s${RESET}\n", helpCommand, helpMessage; \
		} \
	} \
	{ lastLine = $$0 }' $(MAKEFILE_LIST)

goto_app_src:
	-cd ~/labaggregator-bot/

#  Run development environment
dev:
	-BOT_TOKEN=${DEVELOPMENT_BOT_TOKEN} $(DOCKER_COMPOSE_DEV_CMD) up --build $(c)

#  Builds docker compose file in this directory with prod token and launches bot
env_prod: mongod_stop goto_app_src pull
	-BOT_TOKEN=${PRODUCTION_BOT_TOKEN} $(DOCKER_COMPOSE_CMD) up --build -d $(c)

#  Run development environment
env_test: 
	-BOT_TOKEN=${DEVELOPMENT_BOT_TOKEN} $(DOCKER_COMPOSE_TEST_CMD) up --build --abort-on-container-exit -d $(c)

#  Builds docker compose file in this directory with prod token, after git scrip (git pull) command is triggered by github actions
prod_remote: mongod_stop goto_app_src stop
	-${GITHUB_SCRIPT} && BOT_TOKEN=${PRODUCTION_BOT_TOKEN} $(DOCKER_COMPOSE_CMD) up --build -d $(c)

#  stops application
stop: mongod_stop goto_app_src
	-docker compose stop

#  build container
build:
	-$(DOCKER_COMPOSE_CMD) build $(c)

#  up application with git pull
pull:
	-git pull

#  run in daemon mode
up:
	-$(DOCKER_COMPOSE_CMD) up -d $(c)

#  build and run in detached mode + 
buildup:
	-BOT_TOKEN=${DEVELOPMENT_BOT_TOKEN} $(DOCKER_COMPOSE_CMD) up --build 
	
#	&& docker image prune -f $(c)

#  pull from git and reload running containers
pull_n_up: pull up

start:
	-$(DOCKER_COMPOSE_CMD) start $(c)

down:
	-$(DOCKER_COMPOSE_CMD) down $(c)

#  destroy application (containers), preserve images
destroy:
	-$(DOCKER_COMPOSE_CMD) down -v $(c)

#  This will delete all containers, volumers and images
prune: 
	-docker system prune --all --force --volumes 

#  restart application
restart: stop up

#  logs for every container listed in docker-compose.yml
logs:
	-$(DOCKER_COMPOSE_CMD) logs --tail=100 -f $(c)

#  lists application for bot container
logs_bot:
	-$(DOCKER_COMPOSE_CMD) logs --tail=100 -f | grep vendorbot_container $(c)

#  lists containers in docker-compose
ps:
	-$(DOCKER_COMPOSE_CMD) ps

#  gives shell to mongodb_api container
mongodb_shell:
	-docker exec -ti mongodb /bin/bash

#  gives shell to bot container
bot_shell:
	-docker exec -ti vendorbot_container /bin/bash

#  Do update if needed
update_if_needed:
	-git remote update
ifeq ($(shell git rev-parse @), $(shell git rev-parse "origin"))
	-echo "all good"
else
	-echo "Local version is not equalt to remote version"
	-git fetch --all
	-git reset --hard origin/main
	-$(DOCKER_COMPOSE_CMD) up --build -d
endif

#  Install updater to crontab
cron_install_update:
	-crontab -l | grep -q update_if_needed || \
	({ crontab -l; echo "*/15 * * * * cd $(shell pwd -P) && make update_if_needed"; } | crontab -)

#  annihilate all docker images on device
annihilate: 
	docker rm -vf $(shell docker ps -aq) && docker rmi -f $(shell docker images -aq)
