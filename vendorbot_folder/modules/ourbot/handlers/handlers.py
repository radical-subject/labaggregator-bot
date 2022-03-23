import logging
logger = logging.getLogger(__name__)

class Handlers:
    INITIAL = 1
    ADMIN = 1
    LABS = 2
    WISHLIST = 2
    CATEGORIES = 2

    def __init__(self, db_instances):
        self.db_instances=db_instances
        self.db_clients={}
        if db_instances: # If db_instances is not empty.
            for key, value in db_instances.items():
                setattr(self, key, value)
                self.db_clients[f'{key}_client']=value.client
            for key, value in self.db_clients.items():
                setattr(self, key, value)
        self.collection = None
        