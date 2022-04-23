
from typing import Text, List, Tuple
import uuid
import time
from modules.db.blacklist import blacklist_engine
from modules.ourbot.service.cas_to_smiles import is_cas_number
from modules.ourbot.service import batch

import logging
import traceback
logger = logging.getLogger(__name__)


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


def get_contact(user):
    if user["username"]:
        return user["username"]
    elif user["phone_number"]:
        return user["phone_number"]  # TODO добавить + если 7


def reagent_contact(user, reagent):
    if 'contact' in reagent and reagent['contact']:
        return reagent['contact']
    else:
        return get_contact(user)


def find_reagent(user, value):
    reagents = []
    if 'user_reagents' in user:
        for reagent in user['user_reagents']:
            if value in reagent.values():
                reagents.append(reagent)
    return reagents


def get_reagent_contacts(users, text):
    """
    :param users: список объектов пользователей
    :param text: CAS или SMILES
    :return: список контактов
    """
    ret = []
    if users:
        for user in users:
            reagents = find_reagent(user, text)
            for r in reagents:
                contact = reagent_contact(user, r)
                if contact:
                    ret.append(contact)
    return ret


def parse_cas_list(cas_list: List[str], contact: str = ''):
    """
    Фильтруем список CAS, ищем SMILES, удаляем прекурсоры, возвращаем список компонентов для БД и статистику
    :param cas_list:
    :param contact:
    :return:
    """
    valid_cas_list = [r for r in cas_list if is_cas_number(r)]
    failed_cas = [r for r in cas_list if not is_cas_number(r)]
    cas_smiles_list = batch.batch_cas_to_smiles(valid_cas_list)

    no_smiles_list = [cas_smiles[0] for cas_smiles in cas_smiles_list if not cas_smiles[1]]

    cas_smiles_list = [cas_smiles for cas_smiles in cas_smiles_list if cas_smiles[1]]

    cas_smiles_whitelist = []
    errors = []
    for cas, smiles in cas_smiles_list:
        try:
            if not blacklist_engine.is_similar(smiles):
                cas_smiles_whitelist.append((cas, smiles))
        except Exception as err:
            tb = traceback.format_exc()
            logger.error(f"is_similar failed for ({cas}, {smiles}). Error: {tb}")
            errors.append(f"{cas}, {smiles}")

    reagents = []

    now = time.strftime("%d.%m.%Y %H:%M", time.localtime())

    for cas, smiles in cas_smiles_whitelist:
        r = {
            "reagent_internal_id": uuid.uuid4().hex,
            "CAS": cas,
            "SMILES": smiles,
            "sharing_status": "shared",
            "timestamp": now
        }
        if contact:
            r["contact"] = contact
        reagents.append(r)

    message = f"file was successfully parsed and uploaded.\n"
    message += f"<b>import results</b>:\n"
    message += f"Строк в вашем списке <b>{len(cas_list)}</b>\n"
    message += f"Правильных CAS-номеров <b>{len(valid_cas_list)}</b>\n"
    message += f"Опечатка в CAS: <b>{', '.join(failed_cas)}</b>\n"
    message += f"Не найдено SMILES для: <b>{len(no_smiles_list)}</b> позиций\n"
    if no_smiles_list:
        message += "\n".join(no_smiles_list) + "\n"
    message += f"Ошибка обработки SMILES <b>{len(errors)}</b> позиций\n"
    if errors:
        message += "\n".join(errors) + "\n"
    message += f"Найдено SMILES для: <b>{len(cas_smiles_list)}</b> реагентов\n"
    message += f"Прекурсоров найдено и вычеркнуто: <b>{len(cas_smiles_list) - len(cas_smiles_whitelist)}</b>\n"

    return reagents, message


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
