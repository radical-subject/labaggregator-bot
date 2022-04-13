
from typing import Text, List
import uuid
import os, time
os.environ['TZ'] = 'Europe/Moscow'

from modules.ourbot.logger import logger

from modules.ourbot.service.helpers import is_CAS_number

from modules.db.rdkitdb import similarity_search, convert_to_smiles_and_get_additional_data
from modules.ourbot.service.resolver import get_SMILES, batch_SMILES_resolve, CIRPY_resolve


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
                reagent_name: "something 4-something"
                sharing_status: "shared"
            }
        ]
    }

    """
    def __init__(self, *args, **kwargs):
        if args: # If args is not empty.
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

    def create_new_list_reagents(self, CAS_list=list):
        """
        вспомогательная функция, набивающая лист правильно форматированными реагентами для записи в объект
        """
        valid_CAS_numbers = []
        non_valid_CAS_numbers = []

        for CAS_number in CAS_list:
            new_reagent = {
                "reagent_internal_id": uuid.uuid4().hex, 
                "CAS": CAS_number
            }
            if is_CAS_number(CAS_number):
                valid_CAS_numbers.append(new_reagent)
            else:
                non_valid_CAS_numbers.append(CAS_number)

        return valid_CAS_numbers, non_valid_CAS_numbers

    def add_list_of_reagents(self, user_id, contact_username, client, db_instance, CAS_list: List[str]):
        """
        эта функция создает список реагентов в объекте из листа, дополняя его SMILES и потом отсеивая наркотики
        """
        if not CAS_list:
            return

        input_lines_number = len(CAS_list)
        CAS_checker_out = self.create_new_list_reagents(CAS_list)
        
        try:
            self.user_reagents += CAS_checker_out[0]
        except AttributeError:
            setattr(self, "user_reagents", CAS_checker_out[0])

        SMILES_resolver_result = self.resolve_CAS_to_SMILES(contact_username) # дополнение списка SMILES-ами
        blacklist_filter_result = self.blacklist_filter(client, db_instance) # фильтрация структур

        return {
            "input_lines_number": input_lines_number,
            "valid_CAS_numbers": CAS_checker_out[0], 
            "failed_CAS_check_number": CAS_checker_out[1],
            "SMILES_not_found": SMILES_resolver_result["SMILES_not_found"],
            "SMILES_found": SMILES_resolver_result["SMILES_found"],
            "blacklist_filter_result": blacklist_filter_result,
            "total_reagents_imported": len(CAS_checker_out[0])-blacklist_filter_result,
            "total_reagents_count_in_DB": len(self.user_reagents)
        }

    def resolve_CAS_to_SMILES(self, contact_username):
        """
        эта функция читает список реагентов, записи из листа которым не приписаны SMILES, пихает их в отдельный лист.
        кормит этим листом реагентов batch_SMILES_resolve
        полученный лист объектов теперь содержит key "SMILES".
        его дальше дополняет контактными данными и ставит заглушку sharing_status
        """
        
        # забираем в отдельный лист реагенты без CAS
        reagents_without_SMILES_list = [entry for entry in self.user_reagents if not ("SMILES" in entry.keys())]
        # удаляем забранные для парсинга записи
        self.user_reagents = [entry for entry in self.user_reagents if ("SMILES" in entry.keys())]

        logger.info(f'without smiles: {reagents_without_SMILES_list}')
        # парсим записи, добавляя им SMILES
        resolved_list, not_found_list = batch_SMILES_resolve(reagents_without_SMILES_list)

        logger.info(f'resolved_list: {resolved_list}, not_found_list: {not_found_list}')

        # чтобы была возможность в будущем удалить батчем добавленные в одну операцию записи ставим временную метку
        timestamp = time.strftime("%d.%m.%Y %H:%M", time.localtime())
        current_time = timestamp

        # дополняем инфой о контакте, штамп времени и статус шеринга
        for entry in resolved_list:
            entry["contact"] = contact_username
            entry["sharing_status"] = "shared"
            entry["timestamp"] = current_time

        # прибавляем к концу листа реагентов юзера лист с новыми, резолвнутыми реагентами
        self.user_reagents += resolved_list

        return {
            "SMILES_not_found": len(not_found_list),
            "SMILES_found": len(resolved_list)-len(not_found_list)
        }

    def blacklist_filter(self, client, db_instance):
        """
        эта функция удаляет записи из импорта, которые похожи на наркотики
        """
        iterator = self.user_reagents
        drugs_counter = 0
        for entry in iterator:
            SMILES_input = entry["SMILES"]
            SMILES_input = SMILES_input.replace("|", "")
            if SMILES_input != 'resolver_error':
                try:
                    result = similarity_search (client, db_instance, SMILES_input)
                    if result[0][0] > 0.75:
                        # print(f"{convert_to_smiles_and_get_additional_data(client, db_instance, result)[1]['NameRUS']} найден в импорте и вычеркнут")
                        logger.info(f"{entry['CAS']}, SIMILARITY RESULT = {result[0][0]}, {convert_to_smiles_and_get_additional_data(client, db_instance, result)[1]['NameRUS']} найден в импорте и вычеркнут")
                        self.user_reagents.remove(entry)
                        drugs_counter+=1 
                except Exception as e:
                    """
                    вертикальная черта в SMILES - непонятно что несёт, и RDKIT ее не понимает, убираем ее
                    """
                    logger.info(e)
                    pass

        return drugs_counter

    def get_contacts_for_CAS(self, input_CAS):
        contacts = []
        for reagent in self.user_reagents:
            if input_CAS in reagent.values():
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
