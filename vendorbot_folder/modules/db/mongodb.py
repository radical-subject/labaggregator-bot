import os
from datetime import date
import logging
logger = logging.getLogger(__name__)


def dump_database(username, password):
    """
    function for dumping the database. needs root access
    """
    today = date.today()
    # dd.mm.YY
    current_date = today.strftime("%d.%m.%Y")
    path = os.getcwd()
    path = os.path.join(path, "mongodumps", "{}".format(current_date))
    logger.info(path)
    # exec_into_docker_command = "docker exec -it mongodb bash"
    # logging.info(os.system(exec_into_docker_command))

    # запуск команды по дампу
    command = "mongodump --host {} -u {} -p {} --authenticationDatabase admin -o={}"\
        .format("mongodb_api", username, password, path)
    response = os.system(command)
    logger.info(response)

    # запуск команды архивирования
    command = f"zip -r {path}.zip {path}"
    response = os.system(command)
    logger.info(response)

    return response, path


def purge(client, db_instance):
    db_name = db_instance.DATABASE_NAME
    db = client[db_name]
    db.command("dropDatabase")
    logger.info(f"database {db_name} dropped.")


def add_records(client, db_instance, collection_name: str, data: dict):
    db_name = db_instance.DATABASE_NAME
    # logger.info(f"{client}, {db_instance}, {collection_name}, {data}")
    collection = client[db_name][collection_name]
    result = collection.insert_one(data)
    logger.info(f"data inserted into {db_name}, {collection_name}")
    return result


def update_record(client, db_instance, collection_name: str, query: dict, data: dict):
    db_name = db_instance.DATABASE_NAME
    collection = client[db_name][collection_name]
    # result = collection.insert_one(data)
    result = collection.update(query, {"$set": data}, upsert=True)
    logger.info(f"data upserted into {db_name}, {collection_name}")
    return result


def get_records(client, db_instance, collection_name: str, query: dict, *args):
    db_name = db_instance.DATABASE_NAME
    collection = client[db_name][collection_name]
    search = list(collection.find(query, *args))
    # logger.info(search)
    return search
