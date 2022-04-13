import pytest

from modules.ourbot.logger import logger
from modules.db.dbconfig import db_client, MONGO_TEST_DBNAME
from modules.db.dbmodel import UsersCollection


@pytest.fixture
def purge_users_collection() -> None:
    db_client[MONGO_TEST_DBNAME].drop_collection('users_collection')
    logger.info(f"fixture purge: users_collection cleaned.")


@pytest.fixture
def users_collection() -> UsersCollection:
    return UsersCollection(db_client, MONGO_TEST_DBNAME)
