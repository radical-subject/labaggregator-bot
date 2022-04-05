import pymongo, logging


class MongoDriver:
    client = None

    def __init__(self, host, username, password, name):
        try:
            self.DATABASE_HOST = host
            self.DATABASE_ADMIN_USERNAME = username
            self.DATABASE_ADMIN_PASSWORD = password
            logging.info(f"connect to: {self.DATABASE_ADMIN_USERNAME}, {self.DATABASE_ADMIN_PASSWORD}, {self.DATABASE_HOST}")

            self.client = pymongo.MongoClient(self.DATABASE_HOST,
                                              username=self.DATABASE_ADMIN_USERNAME,
                                              password=self.DATABASE_ADMIN_PASSWORD,
                                              authSource="admin")
            # deprecated
            # self.client_base = pymongo.MongoClient(self.DATABASE_HOST)
            # self.client.admin.authenticate(self.DATABASE_ADMIN_USERNAME, self.DATABASE_ADMIN_PASSWORD)
            # self.client = self.client_base(username=self.DATABASE_ADMIN_USERNAME, password=self.DATABASE_ADMIN_PASSWORD)

            self.DATABASE_NAME = name
            logging.info(f"[+] {self.DATABASE_NAME} database connected!")
        except Exception as e:
            logging.info("[-] Database connection error!")
            raise e
