import os
import sys
from logging import getLogger, getLevelName, root, Formatter, StreamHandler

from modules.token_extractor import token
from modules.db.dbconfig import vendorbot_db, blacklist_rdkit_db
from modules.ourbot.ourbot import BotObject

logger = getLogger()


def setup_logger():
    logger.setLevel(getLevelName(os.getenv("LOG_LEVEL", "INFO").upper()))

    formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    console = sys.stderr
    if console is not None:
        console = StreamHandler(console)
        console.setFormatter(formatter)
        logger.addHandler(console)

    logger.info('log level=' + getLevelName(root.level))


def main():

    setup_logger()
    db_instances = dict(
        vendorbot_db=vendorbot_db,
        blacklist_rdkit_db=blacklist_rdkit_db
        )
        
    bot = BotObject(token, **db_instances)
    bot.start()


if __name__ == '__main__':
    main()
