import logging
from modules.db import dbschema
from modules.ourbot.service.decorators import log_errors
logger = logging.getLogger(__name__)


@log_errors
def purge(client, db_instance):
    db_name = db_instance.DATABASE_NAME
    db = client[db_name]
    db.command("dropDatabase")
    logger.info(f"database {db_name} dropped.")

# @log_errors
def add_records(client, db_instance, collection_name: str, data: dict):
    db_name = db_instance.DATABASE_NAME
    # logger.info(f"{client}, {db_instance}, {collection_name}, {data}")
    collection = client[db_name][collection_name]
    result = collection.insert_one(data)
    logger.info(f"data inserted into {db_name}, {collection_name}")
    return result

@log_errors
def update_record(client, db_instance, collection_name: str, query: dict, data: dict):
    db_name = db_instance.DATABASE_NAME
    collection = client[db_name][collection_name]
    # result = collection.insert_one(data)
    result = collection.update(query, {"$set":data}, upsert=True)
    logger.info(f"data upserted into {db_name}, {collection_name}")
    return result

@log_errors
def get_records(client, db_instance, collection_name: str, query: dict, *args):
    db_name = db_instance.DATABASE_NAME
    collection = client[db_name][collection_name]
    search = list(collection.find(query, *args))
    # logger.info(search)
    return search


@log_errors
def get_timerdata_object(client, db_instance, collection_name: str, query: dict, user_id):
        
    # достаем ее из бд
    previous_records=get_records(client, db_instance, collection_name, query)
    
    # результат поиска может оказаться пустым
    if previous_records == [] or previous_records == None:
        timer_object = dbschema.TimerData(
            **{
                "user_id": user_id
            }
        )
        
    else:
        # если раньше у пользователя были записи то импортируем данные пользователя в объект таймера
        timer_object = dbschema.TimerData(
            **previous_records[0]
        )
    
    return timer_object