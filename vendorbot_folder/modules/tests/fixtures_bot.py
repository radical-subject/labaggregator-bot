import pytest
from modules.ourbot.logger import log
from modules.db.dbconfig import vendorbot_db, blacklist_rdkit_db
from modules.ourbot.ourbot import BotObject
from telegram_bot_unittest.routes import TELEGRAM_URL
from telegram_bot_unittest.core import core
from telegram_bot_unittest.user import BOT_TOKEN, Tester, UserBase, ChatBase
from telegram_bot_unittest.fixtures import u
from modules.ourbot.handlers.helpers import LIST_OF_ADMINS


@pytest.fixture(scope='session')
def bot(telegram_server):

    db_instances = dict(
        vendorbot_db=vendorbot_db,
        blacklist_rdkit_db=blacklist_rdkit_db
    )

    log.info('Starting bot for unit-test...')
    bot = BotObject(BOT_TOKEN, TELEGRAM_URL, **db_instances)
    bot.update_dispatcher()
    bot.updater.start_polling()

    yield bot

    bot.updater.stop()


@pytest.fixture(scope='session')
def admin() -> Tester:

    admin = UserBase(LIST_OF_ADMINS[0])
    admin_chat = ChatBase(LIST_OF_ADMINS[0])

    admin = Tester(core, admin, admin_chat)
    return admin
