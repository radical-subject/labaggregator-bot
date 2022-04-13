import pytest

from modules.ourbot.logger import logger
from modules.db.dbconfig import db_client, MONGO_VENDORBOT_DATABASE, MONGO_TEST_DB
from modules.db.dbmodel import UsersCollection


@pytest.fixture
def purge_users_collection() -> None:
    assert MONGO_VENDORBOT_DATABASE == MONGO_TEST_DB, 'set MONGO_VENDORBOT_DATABASE=test_db!'

    db_client[MONGO_VENDORBOT_DATABASE].drop_collection('users_collection')
    logger.info(f"fixture purge: users_collection cleaned.")


@pytest.fixture
def users_collection() -> UsersCollection:
    return UsersCollection(db_client, MONGO_VENDORBOT_DATABASE)
