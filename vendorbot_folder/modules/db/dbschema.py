import re
from typing import Text
import uuid, json 
import os, time, logging
os.environ['TZ'] = 'Europe/Moscow'

from telegram.ext.dispatcher import run_async

from modules.ourbot.service.helpers import is_CAS_number

from modules.db.rdkitdb import similarity_search, convert_to_smiles_and_get_additional_data
from modules.ourbot.service.resolver import get_SMILES, batch_SMILES_resolve, CIRPY_resolve
from modules.ourbot.handlers.decorators import log_errors
logger = logging.getLogger(__name__)


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

    def add_list_of_reagents(self, user_id, contact_username, client, db_instance, CAS_list=list):
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
            "input_lines_number":input_lines_number, 
            "valid_CAS_numbers": CAS_checker_out[0], 
            "failed_CAS_check_number":CAS_checker_out[1], 
            "SMILES_not_found":SMILES_resolver_result["SMILES_not_found"], 
            "SMILES_found":SMILES_resolver_result["SMILES_found"],
            "blacklist_filter_result":blacklist_filter_result,
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
        # парсим записи, добавляя им SMILES
        resolved_SMILES = batch_SMILES_resolve(reagents_without_SMILES_list)

        # чтобы была возможность в будущем удалить батчем добавленные в одну операцию записи ставим временную метку
        timestamp = time.strftime("%d.%m.%Y %H:%M", time.localtime())
        current_time = timestamp

        # дополняем инфой о контакте, штамп времени и статус шеринга
        for entry in resolved_SMILES[0]:
            entry["contact"] = contact_username
            entry["sharing_status"] = "shared"
            entry["timestamp"] = current_time

        # прибавляем к концу листа реагентов юзера лист с новыми, резолвнутыми реагентами
        self.user_reagents += resolved_SMILES[0]

        return {
            "SMILES_not_found":len(resolved_SMILES[1]), 
            "SMILES_found":len(resolved_SMILES[0])-len(resolved_SMILES[1])
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

    def get_user_shared_reagents(self):
        return self.user_reagents

    def get_digest_shared_reagents(self, didgest_part=None):
        """
        takes only shared reagents CAS
        plugs contact info

        output = {
                'CAS1': [contact1, contact2, ...]
                'CAS2': [.... etc, etc]
            }

        """
        if didgest_part: 
            pass
        else:
            digest_part = {}

        for reagent in self.user_reagents: 
            if reagent["CAS"] not in digest_part.keys():
                digest_part[f'{reagent["CAS"]}'] = [reagent["contact"]]
            else: 
                if reagent["contact"] not in digest_part[f'{reagent["CAS"]}']:
                    digest_part[f'{reagent["CAS"]}'].append(reagent["contact"])
        
        return digest_part

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
        UserReagents_export = {**{"_id": self._id}, **dict(self)}
        # json.dumps(UserReagents_export) # exports json string (to use it as python object you should convert it by json.loads())
        return UserReagents_export

























































class wallet:
    """
    This is object for manipulating easily with wallet 
    attributes before sending updated data to database. 

    test_wallet = {
        "wallet_id": "1231241random uuid",
        "money": {
            "total": 1,
            "currency": "RUB"
        },
        "time": {
            "last_change_timestamp": "22.22.2020 14:00:35"
        },
        "ownership": {
            "list_of_owners": ["id1", "id2"]
        },
        "history_of_changes": [
            {
                "owner": "id1",
                "operation": [-2500, "RUB"],
                "timestamp": "22.22.2020 14:00:35"
            }
        ]
    }
    test = wallet(**test_wallet)
    """
    def __init__(self, *args, **kwargs):
        if args: # If args is not empty.
            self.args = args
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)

        # присвоить id если его не было 
        if "wallet_id" not in kwargs.keys():
            self.wallet_id = uuid.uuid4().hex
    
    def __iter__(self):
      for attr, value in self.__dict__.items():
          yield attr, value

    def add_or_retract_money(self, user_id, amount):
        timestamp = time.strftime("%d.%m.%Y %H:%M:%S", time.localtime())
        current_time = timestamp
        try:
            self.money["total"] += amount
            self.history_of_changes.append(
                {
                    "owner": f"{user_id}",
                    "operation": [amount, "RUB"],
                    "timestamp": f"{current_time}"
                }
            )
        except AttributeError:
            setattr(self, "money", {"total": amount, "currency": "RUB"})
            setattr(self, "history_of_changes", [{"owner": f"{user_id}", "operation": [amount, "RUB"], "timestamp": current_time}])

    def get_amount(self):
        return self.money["total"]

    def export(self):
        """
        wallet_export = {
            "name": self.name
            "wallet_id": self.wallet_id,
            "money": self.money,  
            "time": self.time,
            "ownership": self.ownership,
            "history_of_changes": self.history_of_changes
        }
        """      
        wallet_export = {**{"_id": self.wallet_id}, **dict(self)}
        return json.dumps(wallet_export) # exports json string (to use it as python object you should convert it by json.loads())


"""
================================================================================================================================
"""



class TimerData:
    """
    This is object for manipulating easily with timer_data 
    attributes before sending updated data to database. 

    timer_data_structure = {
        _id: "13491823401381",
        user_id: "123123"
        categories: ["deutsch", "chemistry"],
        timerdata:
            [
                {
                    timestamp: "22.22.2020 14:00:35",
                    amount_of_time_in_minutes: 126, 
                    comment: "book reading 300 pages",
                    category_name: "deutsch",
                    archived_status: "active"
                },
                {
                    timestamp: "22.22.2020 14:00:35",
                    amount_of_time_in_minutes: 126
                    comment: "book reading 300 pages"
                    category_name: None,
                    archived_status: "archived"
                }
            ]
    }

    initialization:
    t=TimerData(
        **{
            "smth":smth
        }
    )
    """

    @log_errors
    def __init__(self, *args, **kwargs):
        if args: # If args is not empty.
            self.args = args
        if kwargs:
            for key, value in kwargs.items():
                setattr(self, key, value)

        # присвоить id если его не было 
        if "_id" not in kwargs.keys():
            self._id = uuid.uuid4().hex
    

    @log_errors
    def __iter__(self):
        for attr, value in self.__dict__.items():
            yield attr, value


    @log_errors
    def add_timerdata_entry(self, elapsed_time, comment, category, archived_status: Text, **kwargs):
        timestamp = time.strftime("%Y%m%d%H%M%S", time.localtime())
        current_time = timestamp

        data_dict = {
                    "timerdata_entry_id": uuid.uuid4().hex,
                    "timestamp": f"{current_time}",
                    "amount_of_time_in_minutes": elapsed_time,
                    "comment": comment,
                    "category_name": category,
                    "archived_status": archived_status                  
                }

        if kwargs:       
            data_dict.update(kwargs)

        logger.info(data_dict)

        try:
            self.timerdata.append(data_dict)
        except AttributeError as e:
            logger.info(f"{e}\ncreating...")
            setattr(self, "timerdata", [data_dict])


    @log_errors
    def get_netto_today(self):
        netto_time = 0
        try:
            for i in self.timerdata:
                current_date = time.strftime("%Y%m%d", time.localtime())
                if i["timestamp"][:8] == current_date:
                    netto_time += i["amount_of_time_in_minutes"]
        except AttributeError as e:
            pass
        return netto_time


    @log_errors
    def mutate_archive_status(self, category):
        try:
            if self.timerdata:
                for entry_No in range(len(self.timerdata)):
                    if self.timerdata[entry_No]['category_name'] == category:
                        self.timerdata[entry_No]["archived_status"] = "True" if self.timerdata[entry_No]["archived_status"] == "False" else "False"
                return True
        except:
            return "some error idk"

    @log_errors
    def delete_category(self, category):
        # first remove from lists
        try:
            self.archived_categories.remove(f"{category}")
        except:
            self.categories.remove(f"{category}")
        # then clean up timerdata entries:
        try:
            if self.timerdata:
                for entry_No in range(len(self.timerdata)):
                    if self.timerdata[entry_No]['category_name'] == category:
                        self.timerdata.pop(entry_No)
                return True
        except:
            return "some error idk"

    @log_errors
    def edit_category(self, timerdata_entry_id, category):
        result = next((entry for entry in self.timerdata if entry['timerdata_entry_id'] == timerdata_entry_id), None)
        entry_index = self.timerdata.index(result)
        result['category_name'] = category
        self.timerdata[entry_index] = result

    @log_errors
    def edit_comment(self, timerdata_entry_id, comment):
        result = next((entry for entry in self.timerdata if entry['timerdata_entry_id'] == timerdata_entry_id), None)
        entry_index = self.timerdata.index(result)
        result['comment'] = comment
        self.timerdata[entry_index] = result

    @log_errors
    def edit_timestamp(self, timerdata_entry_id, timestamp):
        result = next((entry for entry in self.timerdata if entry['timerdata_entry_id'] == timerdata_entry_id), None)
        entry_index = self.timerdata.index(result)
        result['timestamp'] = timestamp
        self.timerdata[entry_index] = result

    @log_errors
    def edit_elapsed_time(self, timerdata_entry_id, elapsed_time):
        result = next((entry for entry in self.timerdata if entry['timerdata_entry_id'] == timerdata_entry_id), None)
        entry_index = self.timerdata.index(result)
        result['amount_of_time_in_minutes'] = elapsed_time
        self.timerdata[entry_index] = result

    @log_errors
    def export_json(self):
        """
        timer_data_export = {
        "_id": self._id,
        "user_id": self.user_id,
        "categories": self.categories,
        "archived_categories": self.archived_categories,
        "timerdata": self.timerdata
        }
        """
        timer_data_export = {**dict(self)}
#         timer_data_export = {**{"_id": self._id}, **dict(self)}
        return json.dumps(timer_data_export) # exports json string (to use it as python object you should convert it by json.loads())

    @log_errors
    def export(self):
        """
        timer_data_export = {
        "_id": self._id,
        "user_id": self.user_id,
        "categories": self.categories,
        "archived_categories": self.archived_categories,
        "timerdata": self.timerdata
        }
        """
        timer_data_export = {**dict(self)}
        logger.info("timer_data exported as dict...")
        return timer_data_export # exports dictionary (pymongo uses dict objects as input)


