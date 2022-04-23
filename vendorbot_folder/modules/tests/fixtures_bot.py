import os
import pytest
import logging
from unittest import mock

from modules.bot.bot import BotObject
from modules.bot.helpers import LIST_OF_ADMINS

from modules.tests.routes import TELEGRAM_URL
from modules.tests.core import core
from modules.tests.user import BOT_TOKEN, Tester, UserBase, ChatBase

logger = logging.getLogger(__name__)


@pytest.fixture(autouse=True, scope='session')
def mock_async():
    """
    disable async handlers, cause pytest got error
    :return:
    """
    with mock.patch.dict(os.environ, {"RUN_ASYNC": "false"}):
        yield


@pytest.fixture(scope='session')
def bot(telegram_server,
        mock_async):

    logger.info('Starting Bot for unit-tests...')
    bot = BotObject(BOT_TOKEN, TELEGRAM_URL)
    bot.update_dispatcher()
    bot.updater.start_polling()

    yield bot

    logger.info('Bot stopping...')
    bot.updater.stop()
    logger.info('Bot stopped')


@pytest.fixture(scope='session')
def admin() -> Tester:

    admin = UserBase(LIST_OF_ADMINS[0])
    admin_chat = ChatBase(admin)

    admin = Tester(core, admin, admin_chat)
    return admin


@pytest.fixture(scope='session')
def anonim() -> Tester:
    """
    User without username to test Share contact
    :return: Tester
    """
    u = UserBase()
    u.username = None

    chat = ChatBase(u)

    unknown_user = Tester(core, u, chat)
    return unknown_user
