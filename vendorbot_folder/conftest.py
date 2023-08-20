
import pytest
from modules.tests.typehints import Session


@pytest.hookimpl()
def pytest_sessionstart(session: Session) -> None:

    print("Start pytest testing")


# Load fixtures
pytest_plugins = [
    'telegram_bot_unittest.fixtures',
    'modules.tests.fixtures_db',
    'modules.tests.fixtures_bot',
    'modules.tests.fixtures_cas'
]
