import os
from typing import Dict
from modules.db.dbschema import get_contact
from modules.ourbot.handlers.manage_dialog import get_contact_from_cas_file

FILES_DIR = os.path.join(os.path.dirname(__file__), 'files')


def test_manage(purge_users_collection: None,  # очищаем БД
                bot, user, admin):

    # чтобы создаем пользователя в БД
    user.send_command('/start')
    message = user.get_message()  # убираем из очереди приветственное сообщение

    user.send_command('/manage')

    message = user.get_message()
    assert message
    assert 'Отправьте мне' in message['text']

    # Отправляем файл
    user.send_file(FILES_DIR, 'cas.txt')

    message = user.get_message()
    assert 'Ожидайте: список обрабатывается' in message['text']

    #message = user.get_message(timeout=15.0)
    #assert 'file was successfully parsed and uploaded' in message['text'], message['text']


def test_file_with_contact():
    cas_list = [
        '1-1-1'
    ]

    contact, ret_list = get_contact_from_cas_file(cas_list)
    assert not contact
    assert ret_list == cas_list

    cas_contact_list = [
        'reagents_contact:123',
        '1-1-1'
    ]

    contact, ret_list = get_contact_from_cas_file(cas_contact_list)
    assert contact == '123'
    assert ret_list == cas_list


def test_get_contact(dbuser: Dict):

    dbuser["username"] = "u"
    dbuser["phone_number"] = "1"
    assert get_contact(dbuser) == "u"

    dbuser["username"] = ""
    assert get_contact(dbuser) == "1"

    dbuser["phone_number"] = ""
    assert not get_contact(dbuser)
