import pytest

from modules.db.dbconfig import db_client, MONGO_TEST_DBNAME
from modules.db.dbmodel import UsersCollection


@pytest.fixture
def users_collection() -> UsersCollection:
    return UsersCollection(db_client, MONGO_TEST_DBNAME)
