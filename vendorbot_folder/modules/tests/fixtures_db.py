
from typing import Dict
import pytest
import logging

from modules.db.dbconfig import db_client, MONGO_VENDORBOT_DATABASE, MONGO_TEST_DB
from modules.db.users import UsersCollection
from pymongo.errors import ServerSelectionTimeoutError

logger = logging.getLogger(__name__)


@pytest.fixture
def purge_users_collection() -> None:
    assert MONGO_VENDORBOT_DATABASE == MONGO_TEST_DB, 'set MONGO_VENDORBOT_DATABASE=test_db!'

    try:
        db_client[MONGO_VENDORBOT_DATABASE].drop_collection('users_collection')
        logger.info(f"fixture purge: users_collection cleaned.")
    except ServerSelectionTimeoutError:
        raise Exception('Connect to MongoDB error. Run MongoDB.')


@pytest.fixture
def users_collection() -> UsersCollection:
    return UsersCollection(db_client, MONGO_VENDORBOT_DATABASE)


@pytest.fixture
def dbuser(user,
           users_collection: UsersCollection) -> Dict:

    dbuser = users_collection.get_user(user.id)
    if not dbuser:
        userdata = {
            "_id": user.id,
            "user_id": user.id,
            "username": "@username",
            "firstname": "first_name",
            "lastname": "last_name",
            "phone_number": ""
        }
        users_collection.add_user(userdata)
        dbuser = users_collection.get_user(user.id)

    return dbuser
