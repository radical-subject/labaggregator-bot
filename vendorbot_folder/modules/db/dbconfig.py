import os
import pymongo

from modules.db.mongodb import MongoDriver

MONGO_URL = os.getenv("MONGO_URL")
MONGO_INITDB_ROOT_USERNAME = os.getenv("MONGO_INITDB_ROOT_USERNAME")
MONGO_INITDB_ROOT_PASSWORD = os.getenv("MONGO_INITDB_ROOT_PASSWORD")
MONGO_VENDORBOT_DATABASE = os.getenv("MONGO_VENDORBOT_DATABASE")

MONGO_TEST_DB = 'test_db'

MONGO_BOT_USERNAME = os.getenv("MONGO_BOT_USERNAME")
MONGO_BOT_PASSWORD = os.getenv("MONGO_BOT_PASSWORD")

# Клиент к MongoDB
db_client = pymongo.MongoClient(MONGO_URL,
                                username=MONGO_BOT_USERNAME,
                                password=MONGO_BOT_PASSWORD,
                                authSource="admin",
                                connectTimeoutMS=5)


vendorbot_db = MongoDriver(MONGO_URL, MONGO_BOT_USERNAME, MONGO_BOT_PASSWORD, MONGO_VENDORBOT_DATABASE)

#vendorbot_collections = ("users_collection", "vendors_collection", "crude_vendors_data", "laboratories")


# don't use
#rdkit_db = MongoDriver(MONGO_URL, MONGO_INITDB_ROOT_USERNAME, MONGO_INITDB_ROOT_PASSWORD, "rdkit_db")
#rdkit_collections = ("molecules", "mfp_counts", "permutations")


blacklist_rdkit_db = MongoDriver(MONGO_URL, MONGO_INITDB_ROOT_USERNAME, MONGO_INITDB_ROOT_PASSWORD, "blacklist_rdkit_db")

#blacklist_rdkit_collections = ("molecules", "mfp_counts", "permutations")
