import pymongo, logging


class MongoDriver:
    client = None

    def __init__(self, db_dict: dict):
        try:
            self.DATABASE_HOST = db_dict.get('DATABASE_HOST')
            self.DATABASE_ADMIN_USERNAME = db_dict.get('DATABASE_ADMIN_USERNAME')
            self.DATABASE_ADMIN_PASSWORD = db_dict.get('DATABASE_ADMIN_PASSWORD')
            logging.info(f"connecting with following credentials: {self.DATABASE_ADMIN_USERNAME}, {self.DATABASE_ADMIN_PASSWORD}, {self.DATABASE_HOST}")
            # deprecated:
            # self.client = pymongo.MongoClient(self.DATABASE_HOST)
            
            self.client = pymongo.MongoClient(self.DATABASE_HOST, username=self.DATABASE_ADMIN_USERNAME, password=self.DATABASE_ADMIN_PASSWORD)
            
            # deprecated
            # self.client_base = pymongo.MongoClient(self.DATABASE_HOST)
            # self.client.admin.authenticate(self.DATABASE_ADMIN_USERNAME, self.DATABASE_ADMIN_PASSWORD)
            # self.client = self.client_base(username=self.DATABASE_ADMIN_USERNAME, password=self.DATABASE_ADMIN_PASSWORD)
            self.DATABASE_NAME = db_dict.get('DATABASE_NAME')
            logging.info("[+] Database connected!")
        except Exception as e:
            logging.info("[-] Database connection error!")
            raise e
