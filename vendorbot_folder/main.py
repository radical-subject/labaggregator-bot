import logging
from modules.token_extractor import token #extracts bot token
# database credentials extractor import
from modules.db.dbconfig import root_credentials, timerbot_credentials, vendorbot_credentials, blacklist_rdkit_credentials, rdkit_credentials
from modules.db.mongodb import MongoDriver
# bot object
from modules.ourbot.ourbot import BotObject

# Enable logging
logging.basicConfig(
    # filename='my_runtime_log.log', # saving log to filename
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO  # DEBUG
)
logging.info('logger started')
logger = logging.getLogger(__name__)

def main():

    db_instances = dict(
        root=MongoDriver(root_credentials), 
        timerbot_db=MongoDriver(timerbot_credentials), 
        vendorbot_db=MongoDriver(vendorbot_credentials), 
        blacklist_rdkit_db=MongoDriver(blacklist_rdkit_credentials)
        )
        
    bot = BotObject(token, **db_instances)
    bot.start()

if __name__ == '__main__':
    main()
