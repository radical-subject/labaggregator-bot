
from modules.db.dbconfig import db_client, MONGO_VENDORBOT_DATABASE
import logging
logger = logging.getLogger(__name__)


class UsersCollection:

    def __init__(self, client, db):
        self.name = 'users_collection'
        self.collection = client[db][self.name]
        self.client = client

    def get_user(self, user_id: int):
        return self.collection.find_one({"user_id": user_id})

    def add_user(self, data):
        return self.collection.insert_one(data)

    def get_reagents(self, user_id: int):
        user = self.get_user(user_id)
        if "user_reagents" in user:
            return user["user_reagents"]
        return []

    def get_all_users(self):
        return list(self.collection.find({}))

    def update_user(self, user_id: int, user_data):
        logger.info(f"user {user_id} update data")
        if '_id' in user_data:
            del user_data['_id']  # Performing an update on the path '_id' would modify the immutable field '_id'
        result = self.collection.update_one({"user_id": user_id}, {"$set": user_data}, upsert=True)
        logger.info(f"user {user_id} updated")
        return result

    def get_users_by_cas(self, cas: str):
        return self.collection.find({"user_reagents": {'$elemMatch': {'CAS': cas}}})

    def get_users_by_smiles(self, smiles: str):
        return self.collection.find({"user_reagents": {'$elemMatch': {'SMILES': smiles}}})


users_collection = UsersCollection(db_client, MONGO_VENDORBOT_DATABASE)


def purge(client, db_instance):
    db_name = db_instance.DATABASE_NAME
    db = client[db_name]
    db.command("dropDatabase")
    logger.info(f"database {db_name} dropped.")


def add_records(client, db_instance, collection_name: str, data: dict):
    db_name = db_instance.DATABASE_NAME
    # logger.info(f"{client}, {db_instance}, {collection_name}, {data}")
    collection = client[db_name][collection_name]
    result = collection.insert_one(data)
    logger.info(f"data inserted into {db_name}, {collection_name}")
    return result


def update_record(client, db_instance, collection_name: str, query: dict, data: dict):
    db_name = db_instance.DATABASE_NAME
    collection = client[db_name][collection_name]
    # result = collection.insert_one(data)
    result = collection.update(query, {"$set": data}, upsert=True)
    logger.info(f"data upserted into {db_name}, {collection_name}")
    return result


def get_records(client, db_instance, collection_name: str, query: dict, *args):
    db_name = db_instance.DATABASE_NAME
    collection = client[db_name][collection_name]
    search = list(collection.find(query, *args))
    # logger.info(search)
    return search
