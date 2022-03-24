import os
MONGO_URL = os.environ["MONGO_URL"]
MONGO_INITDB_ROOT_USERNAME = os.environ["MONGO_INITDB_ROOT_USERNAME"]
MONGO_INITDB_ROOT_PASSWORD = os.environ["MONGO_INITDB_ROOT_PASSWORD"]
INIT_DATABASE_NAME = os.environ["MONGO_INITDB_DATABASE"]
MONGO_VENDORBOT_DATABASE = os.environ["MONGO_VENDORBOT_DATABASE"]

MONGO_BOT_USERNAME=os.environ["MONGO_BOT_USERNAME"]
MONGO_BOT_PASSWORD=os.environ["MONGO_BOT_PASSWORD"]

# --------------------------------------
# --------------------------------------
# ROOT ACCESS
# --------------------------------------
# --------------------------------------

root_credentials = dict(
    DATABASE_NAME = INIT_DATABASE_NAME,
    DATABASE_HOST = MONGO_URL,
    DATABASE_ADMIN_USERNAME = MONGO_INITDB_ROOT_USERNAME,
    DATABASE_ADMIN_PASSWORD = MONGO_INITDB_ROOT_PASSWORD
)

root_collections = ()

# --------------------------------------
# --------------------------------------
# BOT DB's
# --------------------------------------
# --------------------------------------


vendorbot_credentials = dict(
    DATABASE_NAME = MONGO_VENDORBOT_DATABASE,
    DATABASE_HOST = MONGO_URL,
    DATABASE_ADMIN_USERNAME = MONGO_BOT_USERNAME,
    DATABASE_ADMIN_PASSWORD = MONGO_BOT_PASSWORD
)

vendorbot_collections = ("users_collection", "vendors_collection", "crude_vendors_data", "laboratories")

timerbot_credentials = dict(
    DATABASE_NAME = "timerbot_db",
    DATABASE_HOST = MONGO_URL,
    DATABASE_ADMIN_USERNAME = MONGO_BOT_USERNAME,
    DATABASE_ADMIN_PASSWORD = MONGO_BOT_PASSWORD
)

timerbot_collections = ()

# --------------------------------------
# --------------------------------------
# CHEMOINFORMATICS
# --------------------------------------
# --------------------------------------

rdkit_credentials = dict(
    DATABASE_NAME = "rdkit_db",
    DATABASE_HOST = MONGO_URL,
    DATABASE_ADMIN_USERNAME = MONGO_INITDB_ROOT_USERNAME,
    DATABASE_ADMIN_PASSWORD = MONGO_INITDB_ROOT_PASSWORD
)

rdkit_collections = ("molecules", "mfp_counts", "permutations")

blacklist_rdkit_credentials = dict(
    DATABASE_NAME = "blacklist_rdkit_db",
    DATABASE_HOST = MONGO_URL,
    DATABASE_ADMIN_USERNAME = MONGO_INITDB_ROOT_USERNAME,
    DATABASE_ADMIN_PASSWORD = MONGO_INITDB_ROOT_PASSWORD
)

blacklist_rdkit_collections = ("molecules", "mfp_counts", "permutations")