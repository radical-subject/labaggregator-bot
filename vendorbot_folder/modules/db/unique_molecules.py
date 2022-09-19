
from modules.db.dbconfig import db_client, MONGO_VENDORBOT_DATABASE
import logging
logger = logging.getLogger(__name__)


class UniqueMolecules:

    def __init__(self, client, db):
        self.name = 'unique_molecules_collection'
        self.collection = client[db][self.name]
        self.client = client

    def get_molecule(self, index: str):
        return self.collection.find_one({"index": index})

    def add_molecule(self, data):
        result = self.collection.insert_one(data)
        #assert result.modified_count == 1

    # def get_reagents(self, user_id: int):
    #     user = self.get_user(user_id)
    #     if "user_reagents" in user:
    #         return user["user_reagents"]
    #     return []

    # def get_all_users(self):
    #     return list(self.collection.find({}))

    def update_molecule(self, index: int, moldoc):
        logger.info(f"molecule entry is being updated")
        query = {'index': index}
        # if '_id' in user_data:
        #     del user_data['_id']  # Performing an update on the path '_id' would modify the immutable field '_id'
        insertion_result = self.collection.update_one(query, {"$set": moldoc}, upsert=True) 
        logger.info(f"molecule entry updated: {insertion_result.acknowledged}")
        #return result
        #assert result.modified_count == 1

    # def get_users_by_cas(self, cas: str):
    #     return list(self.collection.find({"user_reagents": {'$elemMatch': {'CAS': cas}}}))

    # def get_users_by_smiles(self, smiles: str):
    #     return list(self.collection.find({"user_reagents": {'$elemMatch': {'SMILES': smiles}}}))


unique_molecules_collection = UniqueMolecules(db_client, MONGO_VENDORBOT_DATABASE)
