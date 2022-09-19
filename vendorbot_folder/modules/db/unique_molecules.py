
from typing import List
import uuid
import time
from modules.db.blacklist import blacklist_engine
from modules.chem.cas_to_smiles import is_cas_number
from modules.chem import batch

import logging
import traceback
logger = logging.getLogger(__name__)



class UniqueReagentsList:
    """
    This is object for manipulating easily with user 
    attributes before sending updated data to database. 

    test_record = {
        _id: "980159954",
        user_id: "980159954",
        username: "@None",
        firstname: "Alex",
        lastname: "Fedorov"
        laboratory: [
            {
                laboratory_object
            },
            {
                laboratory_object
            }
        ]
        reagent_requests: [
            {
                requested_CAS: "50-00-0"
            }
        ],
        user_reagents: [
            {
                CAS: "50-00-0",
                SMILES: "???",
                #reagent_name: "something 4-something"
                sharing_status: "shared",
                contact: "" # если админ добавил
            }
        ]
    }

    """
    def __init__(self, *args, **kwargs):
        if args:
            self.args = args
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)

        # присвоить id если его не было 
        if "_id" not in kwargs.keys():
            self._id = uuid.uuid4().hex
    
    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value

    def get_contacts_for_reagent(self, value):
        """
        find value in reagents:
        {
            CAS: "75-64-9",
            SMILES: "CC(C)(C)N",
            reagent_name: "something 4-something"
            sharing_status: "shared"
        }
        :param value:
        :return:
        """
        contacts = []
        for reagent in self.user_reagents:
            if value in reagent.values():
                if reagent["contact"] not in contacts:
                    contacts.append(reagent["contact"])
        return contacts

    def export(self):
        """
        UserReagents_export = {
            "_id": self._id
            "user_id": self.user_id,
            "username": self.username,  
            "time": self.time,
            "firstname": self.firstname,
            "lastname": self.lastname,
            "laboratory": self.laboratory,
            "reagent_requests": self.reagent_requests,
            "user_reagents": self.user_reagents
        }
        """
        return {**{"_id": self._id}, **dict(self)}
        # json.dumps(UserReagents_export) # exports json string (to use it as python object you should convert it by json.loads())
