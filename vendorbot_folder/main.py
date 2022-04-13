
from modules.token_extractor import token

from modules.db.dbconfig import vendorbot_db, blacklist_rdkit_db

from modules.ourbot.ourbot import BotObject


def main():
    db_instances = dict(
        vendorbot_db=vendorbot_db,
        blacklist_rdkit_db=blacklist_rdkit_db
        )
        
    bot = BotObject(token, **db_instances)
    bot.start()


if __name__ == '__main__':
    main()
