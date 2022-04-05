import logging
from modules.token_extractor import token #extracts bot token

from modules.db.dbconfig import root_db, timerbot_db, vendorbot_db, blacklist_rdkit_db

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
        root=root_db,
        timerbot_db=timerbot_db,
        vendorbot_db=vendorbot_db,
        blacklist_rdkit_db=blacklist_rdkit_db
        )
        
    bot = BotObject(token, **db_instances)
    bot.start()


if __name__ == '__main__':
    main()
