
import pytest
from modules.db.dbmodel import UsersCollection


def test_add_user(purge_users_collection: None,
                  users_collection: UsersCollection):

    assert not users_collection.get_all_users()

    assert not users_collection.get_user(1)

    userdata = {
        "_id": 1,
        "user_id": 1,
        "username": '@username',
        "firstname": 'first_name',
        "lastname": 'last_name'
    }

    users_collection.add_user(userdata)

    user = users_collection.get_user(1)
    assert user
    assert user['_id'] == 1
    assert user['user_id'] == 1
    assert user['username'] == '@username'
    assert user['firstname'] == 'first_name'
    assert user['lastname'] == 'last_name'

    users = users_collection.get_all_users()
    assert len(users) == 1

    userdata = {
        "_id": 2,
        "user_id": 2,
        "username": '@username2',
        "firstname": 'first_name2',
        "lastname": 'last_name2'
    }
    users_collection.add_user(userdata)

    users = users_collection.get_all_users()
    assert len(users) == 2


