
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

    def get_user_by_reagent_inchi_key(self, inchi_key: str):
        '''
        ищет в коллекции пользователей по листу реагентов совпадения уникального id регагента, и возвращает результат
        '''
        query = {
            "user_reagents": { '$elemMatch': { "inchikey_standard": inchi_key }}
        }
        return self.collection.find(query)

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

    # get location for my reagent
    def get_user_id_by_username(self, update):
        return self.collection.find_one({"username": update.message.from_user.username})['user_id']

    def get_location_by_user_and_cas(self, username: str, cas):
        user_id = self.get_user_id_by_username(username)
        result = self.collection.find_one({"user_id": user_id})['user_reagents']
        locations = []
        for each in result:
            if each['CAS'] in cas:
                if 'location' in each.keys():
                    locations.append(each['location'])

        return '\n'.join(set(locations))


users_collection = UsersCollection(db_client, MONGO_VENDORBOT_DATABASE)
