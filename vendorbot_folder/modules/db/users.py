
import logging
import re
from typing import Optional, List, Tuple, Union

from modules.db.dbconfig import db_client, MONGO_VENDORBOT_DATABASE
from modules.db import dbschema
from modules.reagent import Reagent

logger = logging.getLogger(__name__)


class UsersCollection:
    """
    test_record = {
        _id: "980159954",
        user_id: "980159954",
        username: "@None",
        firstname: "Alex",
        lastname: "Fedorov"
        laboratory: [ {
                laboratory_object
            }, {
                laboratory_object
            }
        ]
        reagent_requests: [ {
                requested_CAS: "50-00-0"
            }
        ],
        user_reagents: [ {
                CAS: "50-00-0",
                SMILES: "???",
                sharing_status: "shared",
                contact: "" # если админ добавил
            }, { ... }
        ]
    }
    """
    def __init__(self, client, db):
        self.name = 'users_collection'
        self.collection = client[db][self.name]
        self.client = client

    def get_user(self, user_id: int):
        return self.collection.find_one({"user_id": user_id})

    def add_user(self, data):
        result = self.collection.insert_one(data)
        if not result.acknowledged:
            logger.error(f"add_user error: {result}")

    def get_reagents(self, user_id: int):
        user = self.get_user(user_id)
        return user.get("user_reagents", [])

    def clear_reagents(self, user_id: int):
        user = self.get_user(user_id)
        user["user_reagents"] = []
        self.update_user(user_id, user)

    def add_reagents(self, user_id: int, reagents):
        # TODO! не знаю, что делает "$set": user_data
        user = self.get_user(user_id)
        arr = [r.to_dict() for r in reagents]
        user["user_reagents"] = user.get("user_reagents", []) + arr
        self.update_user(user_id, user)

    def reagents_count(self, user_id: int):
        return len(self.get_reagents(user_id))

    def get_all_users(self):
        return list(self.collection.find({}))

    def update_user(self, user_id: int, user_data):
        logger.info(f"user {user_id} update data")
        if '_id' in user_data:
            del user_data['_id']  # Performing an update on the path '_id' would modify the immutable field '_id'
        result = self.collection.update_one({"user_id": user_id}, {"$set": user_data}, upsert=True)
        if result.modified_count != 1:
            logger.error(f"user {user_id} updated error: {result}")
        else:
            logger.info(f"user {user_id} updated")

    def get_users_by_reagent_field(self, name: str, value: str):
        return list(self.collection.find({"user_reagents": {'$elemMatch': {name: value}}}))

    def get_users_by_cas(self, cas: str):
        return self.get_users_by_reagent_field('CAS', cas)

    def get_users_by_smiles(self, smiles: str):
        return self.get_users_by_reagent_field('SMILES', smiles)

    def get_user_by_inchi_key(self, inchi_key: str):
        return self.get_users_by_reagent_field('inchikey_standard', inchi_key)

    def get_location_by_user_id_and_cas(self, update, cas):
        user_id = update.message.from_user.id
        result = self.collection.find_one({"user_id": user_id})['user_reagents']
        locations = []
        for each in result:
            if each['CAS'] in cas:
                if 'location' in each.keys():
                    locations.append(each['location'])
        return '\n'.join(set(locations))

    def get_reagents_by_field(self, field, value_list: Union[str, List[str]]) -> List[Reagent]:
        if isinstance(value_list, str):  # можно и список и 1
            value_list = [value_list, ]
        reagents = []
        for value in value_list:
            users = self.get_users_by_reagent_field(field, value)
            for user in users:
                reagents.extend(dbschema.find_reagent(user, value))
        return reagents

    def get_reagents_by_cas(self, cas_list: Union[str, List[str]]) -> List[Reagent]:
        return self.get_reagents_by_field('CAS', cas_list)

    def get_reagents_by_smiles(self, smiles_list: Union[str, List[str]]) -> List[Reagent]:
        return self.get_reagents_by_field('SMILES', smiles_list)

    def get_reagents_by_inchi(self, inchikey_list: Union[str, List[str]]) -> List[Reagent]:
        return self.get_reagents_by_field('inchikey_standard', inchikey_list)

    def get_reagents_by_name(self, name_list: Union[str, List[str]]) -> List[Reagent]:
        """
        В коллекции пользователей ищет по листу реагентов вхождение подстроки text_name в значение поля "name",
        и возвращает контыкты всех положительных результатов.
        """
        if isinstance(name_list, str):  # можно и список и 1
            name_list = [name_list, ]

        reagents = []
        for name in name_list:
            query = [
                {"$match": {
                    "user_reagents": {
                        '$elemMatch': {"name": {"$regex": name, '$options': 'xi'}}}}},
                {"$project": {
                    "user_id": 1,
                    "username": 1,
                    "user_reagents": {
                        "$filter": {
                            "input": "$user_reagents",
                            "cond": {
                                "$regexMatch": {
                                    "input": "$$this.name",
                                    "regex": name,
                                    "options": "xi"
                                }}}}}}]

            for user in list(self.collection.aggregate(query)):
                for reagent in user["user_reagents"]:
                    r = Reagent()
                    r.user_id = user["user_id"]
                    r.from_dict(reagent)
                    reagents.append(r)

        return reagents


users_collection = UsersCollection(db_client, MONGO_VENDORBOT_DATABASE)
