## Setup

Установить docker

## Launch 

### Development: Windows & MacOs & Linux

```
docker-compose -f docker-compose.yml -f docker-compose.development.yml up --build
```
or

```
make env_dev
```

### Production: Windows & MacOs & Linux

```
docker-compose -f up --build
```

or 

```
make env_prod
```

### Файл окружения

.env файл в корне

Содержание:

```
BOT_TOKEN=<token>
```

## Запуск юниттестов

```
docker-compose -f docker-compose.yml -f docker-compose.pytest.yml up --build --abort-on-container-exit
```

or 

```
make env_test
```

## Сервер

команда для очистки неиспользуемых контейнеров:
```
docker system prune -f
```
кажется после каждой смены conda env.yml собирается новый образ на 3.5Гб



#### overwriting env variables for docker-compose 
may be done by prepending KEY=ARG to docker-compose command:

_PRODUCTION_BOT_TOKEN=1342744556:AAEOzLo7yHVLS04YFErYzJNh4-UmPhCltCs docker-compose up --build_

if you want to start bot locally with

_secrets/token.secret_

file, you may omit unnesessary _KEY=ARG_, then bot will fall to default behaviour, reading token from secret file.

### Linux
for git cloning of private repo there should be ssh key added. 
- add private key to /.ssh/ via ssh-add ~/.ssh/id_ed1231.
- then contents of public key should be copied to gihub account settings.

make sure to add PUBLIC KEY to authorized_keys on server via  

cat id_ed25519.pub >> ~/.ssh/authorized_keys
#### makefile deployment 
same dockerfile pipeline may be applied, but also Makefile deploy is prepared.
first clone via ssh on remote server:

_git clone git@github.com:radicalsubject/vendor_bot.git_

make command checks if there are any changes in repository, pulls them and rebuilds + starts project:

_make prod_remote GITHUB_SCRIPT='git pull' PRODUCTION_BOT_TOKEN='1342744556:AAEOzLo7yHVLS04YFErYzJNh4-UmPhCltCs'_

stops application:

_make stop_

lists commands acceptable by make:

_make help_

## Project Structure
there are 5 docker images in development.
- vendorbot_container
- mongodb_api
- nginx
- indigo_service
- ketcher_container

Ketcher is microcervice to draw structures for structure based queries. It depends on indigo and nginx.
Also nginx will be needed for ssl/tls in future development.

vendorbot_container is place where all chemoinformatics (except draw-input by ketcher) is done. For now it is done by mongo-rdkit chemoinformatics data cartridge. 

We look forward to migrate to postgresql catridge published by Timur Madzhidov in nearest future: it is because our web application is based on Express+React+Postgresql.

## gitflow
1. берем задачу из трелло. 
2. делаем форк основного репозитория от dev
3. работаем у себя локально, когда доделали задачу - оформляем pull request
4. владелец репозитория мерджит пул реквест в дев. 

## Further Info
- look into DEVNOTES.md for raw details and notes on development process
