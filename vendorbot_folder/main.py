import logging

from modules.ourbot.logger import log
from modules.token_extractor import token #extracts bot token

from modules.db.dbconfig import root_db, timerbot_db, vendorbot_db, blacklist_rdkit_db

from modules.ourbot.ourbot import BotObject


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
