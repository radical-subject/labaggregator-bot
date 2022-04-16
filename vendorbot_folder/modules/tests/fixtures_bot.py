import pytest
from modules.ourbot.logger import log
from modules.db.dbconfig import vendorbot_db, blacklist_rdkit_db
from modules.ourbot.ourbot import BotObject

from modules.ourbot.handlers.helpers import LIST_OF_ADMINS

from modules.tests.routes import TELEGRAM_URL
from modules.tests.core import core
from modules.tests.user import BOT_TOKEN, Tester, UserBase, ChatBase


@pytest.fixture(scope='session')
def bot(telegram_server):

    db_instances = dict(
        vendorbot_db=vendorbot_db,
        blacklist_rdkit_db=blacklist_rdkit_db
    )

    log.info('Starting Bot for unit-tests...')
    bot = BotObject(BOT_TOKEN, TELEGRAM_URL, **db_instances)
    bot.update_dispatcher()
    bot.updater.start_polling()

    yield bot

    log.info('Bot stopping...')
    bot.updater.stop()
    log.info('Bot stopped')


@pytest.fixture(scope='session')
def admin() -> Tester:

    admin = UserBase(LIST_OF_ADMINS[0])
    admin_chat = ChatBase(LIST_OF_ADMINS[0])

    admin = Tester(core, admin, admin_chat)
    return admin
