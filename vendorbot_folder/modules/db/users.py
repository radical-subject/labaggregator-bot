
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
        result = self.collection.insert_one(data)
        #assert result.modified_count == 1

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
        #return result
        #assert result.modified_count == 1

    def get_users_by_cas(self, cas: str):
        return list(self.collection.find({"user_reagents": {'$elemMatch': {'CAS': cas}}}))

    def get_users_by_smiles(self, smiles: str):
        return list(self.collection.find({"user_reagents": {'$elemMatch': {'SMILES': smiles}}}))


users_collection = UsersCollection(db_client, MONGO_VENDORBOT_DATABASE)
