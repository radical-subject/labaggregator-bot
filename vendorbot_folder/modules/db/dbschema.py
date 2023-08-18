
import uuid
import logging

from modules.reagent import REAGENT_SHARED

logger = logging.getLogger(__name__)


def get_shared_reagents(user):
    if 'user_reagents' in user:
        return filter(lambda r: r['sharing_status'] == REAGENT_SHARED, user['user_reagents'])
    return []


def reagent_name(r):
    """
    У нас может быть что-то не заполнено (?)
    """
    if 'CAS' in r and r['CAS']:
        return r['CAS']
    # как добавим, тогда раскомментирую
    # if 'reagent_name' in r and r['reagent_name']:
    #    return r['reagent_name']


def reagent_CAS(r):
    """
    У нас может быть что-то не заполнено (?)
    """
    if 'CAS' in r and r['CAS']:
        return r['CAS']


def get_contact(user):
    if user["username"]:
        return f"@{user['username']}"
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



class UserReagents:
    """
    TODO: удалить этот класс

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
