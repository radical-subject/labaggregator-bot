
import pytest
from typing import Dict
from modules.db.dbmodel import UsersCollection


def test_add_user(purge_users_collection: None,
                  users_collection: UsersCollection):

    assert not users_collection.get_all_users()

    assert not users_collection.get_user(1)

    userdata = {
        "_id": 1,
        "user_id": 1,
        "username": "@username",
        "firstname": "first_name",
        "lastname": "last_name",
        "phone_number": ""
    }

    users_collection.add_user(userdata)

    user = users_collection.get_user(1)
    assert user
    assert user["_id"] == 1
    assert user["user_id"] == 1
    assert user["username"] == "@username"
    assert user["firstname"] == "first_name"
    assert user["lastname"] == "last_name"
    assert user["phone_number"] == ""

    users = users_collection.get_all_users()
    assert len(users) == 1

    userdata = {
        "_id": 2,
        "user_id": 2,
        "username": "@username2",
        "firstname": "first_name2",
        "lastname": "last_name2",
        "phone_number": ""
    }
    users_collection.add_user(userdata)

    users = users_collection.get_all_users()
    assert len(users) == 2


def test_update_user(purge_users_collection: None,
                     users_collection: UsersCollection,
                     dbuser: Dict,
                     user_id: int):

    dbuser["user_reagents"] = [
        {
            "CAS": "1-1-1"
        },
        {
            "CAS": "2-2-2"
        }
    ]
    users_collection.update_user(user_id, dbuser)
    dbuser = users_collection.get_user(user_id)

    assert len(users_collection.get_reagents(user_id)) == 2

    dbuser["phone_number"] = "1"
    users_collection.update_user(user_id, dbuser)

    dbuser = users_collection.get_user(user_id)
    assert dbuser["phone_number"] == "1"

    assert len(users_collection.get_reagents(user_id)) == 2

    dbuser["user_reagents"] = [
        {
            "CAS": "1-1-1"
        },
        {
            "CAS": "2-2-2"
        },
        {
            "CAS": "3-3-3"
        }
    ]

    users_collection.update_user(user_id, dbuser)
    assert len(users_collection.get_reagents(user_id)) == 3
