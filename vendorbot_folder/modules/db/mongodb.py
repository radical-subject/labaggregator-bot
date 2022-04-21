import os
from datetime import date
import pymongo
import logging
logger = logging.getLogger(__name__)


class MongoDriver:
    client = None

    def __init__(self, host, username, password, name):
        try:
            self.DATABASE_HOST = host
            self.DATABASE_ADMIN_USERNAME = username
            self.DATABASE_ADMIN_PASSWORD = password
            logger.info(f"connect to: {self.DATABASE_ADMIN_USERNAME}, {self.DATABASE_ADMIN_PASSWORD}, {self.DATABASE_HOST}")

            self.client = pymongo.MongoClient(self.DATABASE_HOST,
                                              username=self.DATABASE_ADMIN_USERNAME,
                                              password=self.DATABASE_ADMIN_PASSWORD,
                                              authSource="admin")
            # deprecated
            # self.client_base = pymongo.MongoClient(self.DATABASE_HOST)
            # self.client.admin.authenticate(self.DATABASE_ADMIN_USERNAME, self.DATABASE_ADMIN_PASSWORD)
            # self.client = self.client_base(username=self.DATABASE_ADMIN_USERNAME, password=self.DATABASE_ADMIN_PASSWORD)

            self.DATABASE_NAME = name
            logger.info(f"[+] {self.DATABASE_NAME} database connected!")
        except Exception as e:
            logger.info("[-] Database connection error!")
            raise e


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
