import os
import pymongo


MONGO_URL = os.getenv("MONGO_URL")
MONGO_INITDB_ROOT_USERNAME = os.getenv("MONGO_INITDB_ROOT_USERNAME")
MONGO_INITDB_ROOT_PASSWORD = os.getenv("MONGO_INITDB_ROOT_PASSWORD")
MONGO_VENDORBOT_DATABASE = os.getenv("MONGO_VENDORBOT_DATABASE")
MOLECULES_DATABASE = 'molecules_db'

MONGO_TEST_DB = 'test_db'

MONGO_BOT_USERNAME = os.getenv("MONGO_BOT_USERNAME")
MONGO_BOT_PASSWORD = os.getenv("MONGO_BOT_PASSWORD")

# Клиент к MongoDB
db_client = pymongo.MongoClient(MONGO_URL,
                                username=MONGO_BOT_USERNAME,
                                password=MONGO_BOT_PASSWORD,
                                authSource="admin",
                                connectTimeoutMS=5)


root_client = pymongo.MongoClient(MONGO_URL,
                                  username=MONGO_INITDB_ROOT_USERNAME,
                                  password=MONGO_INITDB_ROOT_PASSWORD,
                                  authSource="admin",
                                  connectTimeoutMS=5)

