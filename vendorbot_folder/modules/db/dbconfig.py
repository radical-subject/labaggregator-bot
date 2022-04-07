import os
import pymongo

from modules.db.mongodb import MongoDriver

MONGO_URL = os.environ["MONGO_URL"]
MONGO_INITDB_ROOT_USERNAME = os.environ["MONGO_INITDB_ROOT_USERNAME"]
MONGO_INITDB_ROOT_PASSWORD = os.environ["MONGO_INITDB_ROOT_PASSWORD"]
INIT_DATABASE_NAME = os.environ["MONGO_INITDB_DATABASE"]
MONGO_VENDORBOT_DATABASE = os.environ["MONGO_VENDORBOT_DATABASE"]

MONGO_BOT_USERNAME = os.environ["MONGO_BOT_USERNAME"]
MONGO_BOT_PASSWORD = os.environ["MONGO_BOT_PASSWORD"]

MONGO_TEST_DBNAME = 'test_db'

# Клиент к MongoDB
db_client = pymongo.MongoClient(MONGO_URL,
                                username=MONGO_BOT_USERNAME,
                                password=MONGO_BOT_PASSWORD,
                                authSource="admin")


vendorbot_db = MongoDriver(MONGO_URL, MONGO_BOT_USERNAME, MONGO_BOT_PASSWORD, MONGO_VENDORBOT_DATABASE)

#vendorbot_collections = ("users_collection", "vendors_collection", "crude_vendors_data", "laboratories")


# don't use
#rdkit_db = MongoDriver(MONGO_URL, MONGO_INITDB_ROOT_USERNAME, MONGO_INITDB_ROOT_PASSWORD, "rdkit_db")
#rdkit_collections = ("molecules", "mfp_counts", "permutations")


blacklist_rdkit_db = MongoDriver(MONGO_URL, MONGO_INITDB_ROOT_USERNAME, MONGO_INITDB_ROOT_PASSWORD, "blacklist_rdkit_db")

#blacklist_rdkit_collections = ("molecules", "mfp_counts", "permutations")
