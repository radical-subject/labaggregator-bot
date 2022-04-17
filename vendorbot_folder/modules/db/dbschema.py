
from typing import Text, List, Tuple
import uuid
import time
from modules.db.blacklist import blacklist_engine
from modules.ourbot.service.helpers import is_cas_number
from modules.ourbot.service.cas_to_smiles import banch_cas_to_smiles


def get_shared_reagents(user):
    if 'user_reagents' in user:
        return filter(lambda r: r['sharing_status'] == 'shared', user['user_reagents'])
    return []


def reagent_name(r):
    """
    У нас может быть что-то не заполнено (?)
    """
    if 'CAS' in r and r['CAS']:
        return r['CAS']
    if 'reagent_name' in r and r['reagent_name']:
        return r['reagent_name']


def reagent_CAS(r):
    """
    У нас может быть что-то не заполнено (?)
    """
    if 'CAS' in r and r['CAS']:
        return r['CAS']


def reagent_contact(user, reagent):
    if 'contact' in reagent and reagent['contact']:
        return reagent['contact']
    elif 'username' in user and user['username']:
        return user['username']
    return user['user_id']   # хоть так


def parse_cas_list(cas_list: List[str], contact: str = ''):
    """
    Фильтруем список CAS, ищем SMILES, удаляем прекурсоры, возвращаем список компонентов для БД и статистику
    :param cas_list:
    :param contact:
    :return:
    """
    valid_cas_list = [r for r in cas_list if is_cas_number(r)]
    failed_cas = [r for r in cas_list if not is_cas_number(r)]
    cas_smiles_list = banch_cas_to_smiles(valid_cas_list)

    cas_smiles_whitelist = [cas_smile for cas_smile in cas_smiles_list if
                            not blacklist_engine.is_similar(cas_smile[1])]

    reagents = []

    now = time.strftime("%d.%m.%Y %H:%M", time.localtime())

    for cas, smiles in cas_smiles_whitelist:
        reagents.append({
            "reagent_internal_id": uuid.uuid4().hex,
            "CAS": cas,
            "SMILES": smiles,
            "contact": contact,
            "sharing_status": "shared",
            "timestamp": now
        })

    return reagents, f"""file was successfully parsed and uploaded.
<b>import results</b>:
Строк в вашем списке: <b>{len(cas_list)}</b>
Правильных CAS-номеров: <b>{len(valid_cas_list)}</b>
Опечатка в CAS: <b>{", ".join(failed_cas)}</b>
Не найдено SMILES для: <b>{len(valid_cas_list) - len(cas_smiles_list)}</b> позиций
Найдено SMILES для: <b>{len(cas_smiles_list)}</b> реагентов
Прекурсоров найдено и вычеркнуто: <b>{len(cas_smiles_list) - len(cas_smiles_whitelist)}</b>
"""


class UserReagents:
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
                reagent_name: "something 4-something"
                sharing_status: "shared"
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

    def add_phone_number (self, phone_number):
        self.phone_number = phone_number

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
