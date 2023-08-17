
import logging
import re

from modules.db.dbconfig import db_client, MONGO_VENDORBOT_DATABASE
from modules.db import dbschema

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

    def get_location_by_user_and_cas(self, update, cas):
        user_id = update.message.from_user.id
        result = self.collection.find_one({"user_id": user_id})['user_reagents']
        locations = []
        for each in result:
            if each['CAS'] in cas:
                if 'location' in each.keys():
                    locations.append(each['location'])

        return '\n'.join(set(locations))

    def get_location_by_user_and_inchi_key(self, update, inchi_key: str):
        '''
        ищет в коллекции пользователей нужного пользователя и по листу реагентов - ищет совпадения уникального id регагента.
        возвращает поле локации для данного реагента
        '''
        
        user_id = update.message.from_user.id
        self.get_user_by_reagent_inchi_key(inchi_key)
        query = {
            "user_id": user_id,
            "user_reagents": { '$elemMatch': { "inchikey_standard": inchi_key }}
        }
        results = self.collection.find(query)
        
        locations = []
        for result in results:
            for each in result['user_reagents']:
                if "inchikey_standard" in each.keys():
                    if each["inchikey_standard"] == inchi_key:
                        if 'location' in each.keys():
                            locations.append(each['location'])

        return '\n'.join(set(locations))
    
    def get_reagents_contacts_by_inchi_key(self, inchi_key: str):
        '''
        ищет в коллекции пользователей по листу реагентов совпадения уникального id регагента, и возвращает результат
        '''

        query = {
            "user_reagents": { '$elemMatch': { "inchikey_standard": inchi_key }}
            }
        contacts = []
        users = self.collection.find(query)
        for user in users:
            reagents = dbschema.find_reagent(user, inchi_key)
            
            for reagent in reagents: 
                contacts += [dbschema.reagent_contact(user, reagent)]

        return contacts

    def get_my_reagents_by_text_name(self, text_name, user_id):
        '''
        Ищет у данного пользователя по листу реагентов вхождение подстроки text_name в значение поля "name",
        и возвращает месторасположения (location) всех положительных результатов.
        '''

        query = [
            {"$match": {
                "user_id": user_id, 
                "user_reagents": {
                    '$elemMatch': { "name": { "$regex": text_name, '$options' : 'xi'}}}}},
            {"$project": {
                "user_id" : 1,
                "username": 1,
                "user_reagents": {
                    "$filter": {
                        "input": "$user_reagents",
                        "cond": {
                            "$regexMatch": {
                                "input": "$$this.name",
                                "regex": text_name,
                                "options": "xi"
        }}}}}}]

        my_similar_reagents = list(self.collection.aggregate(query))
        if my_similar_reagents:
            reagents = dict((one['name'], '\n'.join(set(re.split(',|\n', one['location'])))) for one in my_similar_reagents[0]['user_reagents'])
            return reagents

        else:
            return None


    def get_reagents_by_text_name(self, text_name):
        '''
        В коллекции пользователей ищет по листу реагентов вхождение подстроки text_name в значение поля "name",
        и возвращает контыкты всех положительных результатов.
        '''

        query = [
            {"$match": {
                "user_reagents": {
                    '$elemMatch': { "name": { "$regex": text_name, '$options' : 'xi'}}}}},
            {"$project": {
                "user_id" : 1,
                "username": 1,
                "user_reagents": {
                    "$filter": {
                        "input": "$user_reagents",
                        "cond": {
                            "$regexMatch": {
                                "input": "$$this.name",
                                "regex": text_name,
                                "options": "xi"
        }}}}}}]
        
        reagents = dict()
        user_objects_dict = dict((user_object['user_id'], user_object) for user_object in self.collection.aggregate(query))

        for user in user_objects_dict:

            reagents.update({user: ''})

            for reagent in user_objects_dict[user]['user_reagents']:

                new_value = reagents[user] + reagent['name'].replace('\n', ',') + '\n'
                    
                reagents.update({
                   user: '\n'.join(set([name.strip() for name in new_value.replace('\n', ',').split(',')]))
                })

        return reagents


users_collection = UsersCollection(db_client, MONGO_VENDORBOT_DATABASE)
