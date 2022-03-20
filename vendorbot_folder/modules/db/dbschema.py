from typing import Text
import uuid, json 
import os, time, logging
os.environ['TZ'] = 'Europe/Moscow'

from modules.ourbot.service.decorators import log_errors
logger = logging.getLogger(__name__)

time.tzset()

# time.strftime("%d.%m.%Y %H:%M:%S")
# time.strftime("%d.%m.%Y %H:%M:%S", time.localtime())

#
# {
#     _id: "980159954",
#     categories:
#     [
#         {
#             category_name: "chemistry"
#             archived_status: "Active"|"Archived"
#             timer_data:
#             [
#                 {
#                     date: "20200617",
#                     amount_of_time_in_minutes: 126
#                     comment: "book reading 300 pages"
#                 },
#                 {
#                     date: "20200617",
#                     amount_of_time_in_minutes: 126
#                     comment: "book reading 300 pages"
#                 }
#             ]
#         },
#         {
#             category_name: "deutsch"
#             archived_status: "Active"|"Archived"
#             timer_data:
#             [
#                 {
#                     date: "20200617",
#                     amount_of_time_in_minutes: 126
#                     comment: "book reading 300 pages"
#                 },
#                 {
#                     date: "20200617",
#                     amount_of_time_in_minutes: 126
#                     comment: "book reading 300 pages"
#                 }
#             ]
#         }
#     ]
# }

"""

timer_data_structure = {
    _id: "980159954",
    categories: ["deutsch", "chemistry"],
    timer_data:
        [
            {
                date: "20200617",
                amount_of_time_in_minutes: 126, 
                comment: "book reading 300 pages",
                category_name: "deutsch",
                archived_status: "active"
            },
            {
                date: "20200617",
                amount_of_time_in_minutes: 126
                comment: "book reading 300 pages"
                category_name: None,
                archived_status: "archived"
            }
        ]
}

userdata_dict = {
    _id: "980159954",
    user_id: "980159954",
    username: "@None",
    firstname: "Alex",
    lastname: "Fedorov"
}


"""


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


