
from typing import Dict
import pytest
import logging

from modules.db.dbconfig import db_client, MONGO_VENDORBOT_DATABASE, MONGO_TEST_DB
from modules.db.dbmodel import UsersCollection

logger = logging.getLogger(__name__)


@pytest.fixture
def purge_users_collection() -> None:
    assert MONGO_VENDORBOT_DATABASE == MONGO_TEST_DB, 'set MONGO_VENDORBOT_DATABASE=test_db!'

    db_client[MONGO_VENDORBOT_DATABASE].drop_collection('users_collection')
    logger.info(f"fixture purge: users_collection cleaned.")


@pytest.fixture
def users_collection() -> UsersCollection:
    return UsersCollection(db_client, MONGO_VENDORBOT_DATABASE)


@pytest.fixture
def dbuser(user_id: int,
           users_collection: UsersCollection) -> Dict:

    user = users_collection.get_user(user_id)
    if not user:
        userdata = {
            "_id": user_id,
            "user_id": user_id,
            "username": "@username",
            "firstname": "first_name",
            "lastname": "last_name",
            "phone_number": ""
        }
        users_collection.add_user(userdata)
        user = users_collection.get_user(user_id)
        
    return user
