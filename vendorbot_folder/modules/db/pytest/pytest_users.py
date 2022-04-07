
import pytest
from modules.db import dbmodel
from modules.db.dbmodel import UsersCollection


def test_add_user(users_collection: UsersCollection):

    userdata = {
        "_id": 1,
        "user_id": 1,
        "username": '@username',
        "firstname": 'first_name',
        "lastname": 'last_name'
    }

    assert not users_collection.get_user(1)

    users_collection.add_user(userdata)

    user = users_collection.get_user(1)
    assert user
