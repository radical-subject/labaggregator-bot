
import pytest
import logging
logger = logging.getLogger(__name__)
from .routes import start_server, shutdown_server
from .user import Tester, UserBase, ChatBase
from .core import core


u = UserBase()
chat = ChatBase(u)


@pytest.fixture(scope='session')
def user() -> Tester:
    user = Tester(core, u, chat)
    return user


@pytest.fixture(scope='session')
def telegram_server():
    s, t = start_server()
    logger.info('telegram server started')
    yield
    logger.info('telegram server shutdown begin')
    shutdown_server(s, t)
