import os
from typing import Dict
from modules.db.dbschema import get_contact
from modules.ourbot.handlers.manage_dialog import get_contact_from_cas_file
from modules.db.dbmodel import users_collection

FILES_DIR = os.path.join(os.path.dirname(__file__), 'files')


def test_manage(purge_users_collection: None,  # очищаем БД
                bot, user, admin,
                mock_batch_cas_to_smiles,
                mock_is_similar):

    mock_batch_cas_to_smiles.side_effect = [[("2749-11-3", "CC"), ("75-64-9", "CCC"), ("120-46-7", "CCCC")]]
    mock_is_similar.side_effect = [False, False, False]

    # чтобы создаем пользователя в БД
    user.send_command('/start')
    message = user.get_message()  # убираем из очереди приветственное сообщение

    user.send_command('/manage')

    message = user.get_message()
    assert message
    assert 'Отправьте мне' in message['text']

    # Отправляем файл
    user.send_file(FILES_DIR, 'cas1.txt')

    message = user.get_message()
    assert 'Ожидайте: список обрабатывается' in message['text']

    message = user.get_message()

    reagents = users_collection.get_reagents(user.id)
    assert len(reagents) == 3

    text = message['text']
    assert "Строк в вашем списке <b>4</b>" in text
    assert "Правильных CAS-номеров <b>3</b>" in text
    assert "Опечатка в CAS: <b>1-1-1</b>" in text
    assert "Не найдено SMILES для: <b>0</b> позиций" in text
    assert "Ошибка обработки SMILES <b>0</b> позиций" in text
    assert "Найдено SMILES для: <b>3</b> реагентов" in text
    assert "Прекурсоров найдено и вычеркнуто: <b>0</b>" in text
    assert "Итого: База реагентов перезаписана." in text
    assert "Содержит <b>3</b> реагентов" in text

    # после завершения мы показываем приветственное
    message = user.get_message()


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
