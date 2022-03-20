# rewriting of vendorbot aggregator!
proof of concept for:
- fine chemicals for RnD purpose marketplace (aggregator)
- structure and substructure search and chemical drawing
- database blacklisting
- user wishlists and commercial requests by mail

## manual creating environment
* _conda create -n vendorbot python=3.8_
* _conda env export --name vendorbot > env.yml_
or _conda env export --no-builds > env.yml_ - this will skip builds and theoretically make env.yml more crossplatform, but i prefer to write down all dependencies manually 
* _conda env create -f env.yml_
* _conda activate vendorbot
* _conda env list

## manual DB creation and configuration 
#### configure the locally hosted mongodb
windows: _net start/stop MongoDB_
go to mongo shell and type: 

* _use admin_
* _db.createUser({user: "labaggregator", pwd: "", roles: [ { role: "userAdminAnyDatabase", db: "admin" }, "readWriteAnyDatabase" ]})_
* _db.grantRolesToUser( "labaggregator", ["userAdminAnyDatabase"])_
* _db.auth("lolilithium", "")_
* _db.grantRolesToUser( "labaggregator", ["root"])_

## DOCKER STARTUP
docker-compose up --build (don't forget to OPEN Docker App before)

attention:
__token lies in token.secret file which is included in .gitignore with *.secret__

## DOCKER STARTUP
The docker container with bot source named _bot_container_ is configured by Dockerfile.
The docker container with database named _mongodb_api_ is created from docker image: _mongo:latest_

Docker-compose.yml is the file which binds two containers in one app. The bot is run by single command:
_docker-compose build && docker-compose up_
or
_docker-compose up --build_ 
the latter will rebuild the project if any changes were made, but the database data will remain intact. beware!

### for development we create docker-compose.development.yml with some fields that shoul be overriden for dev start: 
_docker-compose -f docker-compose.yml -f docker-compose.development.yml up --build_ 

restarting of the whole project inside docker container is acheved by nodemon, which is triggered by save of any file inside watched folder (bot_container)

also we share volumes with data between mongodb_api and bot_container by following example declaration in docker-compose.yml:

.....
    volumes:
        - ./mongo-init.js:/docker-entrypoint-initdb.d/mongo-init.js:ro
        - mongodb_api:/data/db
    # shorten logs
    command: mongod --quiet # --logpath /dev/null # логи монго выключены насовсем
volumes:
    mongodb_api:
.....

**mongo-init.js** theoretically contains information for prior initialization and addition of new database user, but the login is processed by logging into _admin_ database with 
*MONGO_INITDB_ROOT_USERNAME: root*
*MONGO_INITDB_ROOT_PASSWORD: root* 
specified in _environment variables_ in docker-compose.yml

*MONGO_INITDB_ROOT_USERNAME: root* and *MONGO_INITDB_ROOT_PASSWORD: root* and also **MONGO_URL=mongodb://mongodb_api:27017/bot** are passed as environment variables.
**MONGO_URL=mongodb://mongodb_api:27017/bot** points from bot container to database deployed in database container.
*bot* container depends on *mongodb_api* container and is build and deployed only after database container.

## FOR DUMPING THE DB 
first log into running container, then install mongodb-org-tools, then and then _mongodump_ and _mongorestore_ are available.
_docker exec -it mongodb bash_
_apt update && apt install mongodb-org-tools_

## makefile pipeline
first clone via ssh on remote server 
git clone git@github.com:radicalsubject/vendor_bot.git

there should be ssh key added. private to /.ssh/ via ssh-add ~/.ssh/id_ed1231 
then contents of public key should be copied to gihub account settings 

make sure to add PUBLIC KEY to authorized_keys on server via  
cat id_ed25519.pub >> ~/.ssh/authorized_keys

bot token may be passed in env variable

make prod PRODUCTION_BOT_TOKEN="1342744556:AAEOzLo7yHVLS04YFErYzJNh4-UmPhCltCs"

make stop
stops application

make help - lists commands acceptable by make

## overwriting env variables for docker-compose 
may be done by prepending KEY=ARG to docker-compose command
PRODUCTION_BOT_TOKEN=1342744556:AAEOzLo7yHVLS04YFErYzJNh4-UmPhCltCs docker-compose up --build

if you want to start bot locally with old token.secret file you may omit unnesessary KEY=ARG then bot will fall to default behaviour


