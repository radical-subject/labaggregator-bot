
import time
from modules.db.dbmodel import users_collection
from modules.ourbot.handlers.search_dialog import CANCEL_SEARCH
from modules.ourbot.service.pytests.pytest_cas_to_smiles import PubChempyComponent


reagents1 = [{'reagent_internal_id': '1', 'CAS': '2749-11-3', 'SMILES': 'C[C@H]([NH3+])CO',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'},
             {'reagent_internal_id': '2', 'CAS': '75-64-9', 'SMILES': 'CC(C)(C)N',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'},
             {'reagent_internal_id': '3', 'CAS': '120-46-7', 'SMILES': 'O=C(CC(=O)c1ccccc1)c2ccccc2',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'}]

reagents2 = [{'reagent_internal_id': '4', 'CAS': '106-97-8', 'SMILES': '106-97-8',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'},
             {'reagent_internal_id': '5', 'CAS': '13896-65-6', 'SMILES': '[Ru+3].[I-].[I-].[I-]',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'},
             {'reagent_internal_id': '6', 'CAS': '120-46-7', 'SMILES': 'O=C(CC(=O)c1ccccc1)c2ccccc2',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'}]


def wait_create_user(user_id, timeout=2.0):
    """
    Почему то не сразу появляется пользователь в БД,
    поэтому нужно чуть подождать
    :param user_id:
    :param timeout:
    :return:
    """
    user = None
    start = time.time()
    while not user or time.time() - start < timeout:
        user = users_collection.get_user(user_id)
        time.sleep(0.1)
    return user


def create_user_with_reagents(user, reagents):
    """
    Создаем пользователя в БД, заполняем ему реагенты
    :param user:
    :param reagents:
    :return:
    """
    user.send_command('/start')
    user.get_message()
    dbuser = wait_create_user(user.id)
    dbuser["user_reagents"] = reagents
    users_collection.update_user(user.id, dbuser)
    return dbuser


def test_search_bad_cas(purge_users_collection: None,  # очищаем БД
                        bot, user, admin,
                        mock_pubchempy_get_compounds):

    mock_pubchempy_get_compounds.side_effect = [None, None]

    dbuser = create_user_with_reagents(user, reagents1)
    dbadmin = create_user_with_reagents(admin, reagents2)

    user.send_command("/search")

    text = user.get_message_text()
    assert "Пришли интересующий CAS-номер" in text

    user.send_message("1-1-1")

    text = user.get_message_text()
    assert "Не похоже на CAS. Сейчас поищем по названию..." in text

    text = user.get_message_text()
    assert "Реагентом пока никто не готов поделиться." in text

    user.send_message(CANCEL_SEARCH)
    text = user.get_message_text()
    assert "Поиск завершен" in text


def test_search_cas(purge_users_collection: None,  # очищаем БД
                    bot, user, admin):

    dbuser = create_user_with_reagents(user, reagents1)
    dbadmin = create_user_with_reagents(admin, reagents2)

    user.send_command("/search")

    text = user.get_message_text()
    assert "Пришли интересующий CAS-номер" in text

    user.send_message("2749-11-3")

    text = user.get_message_text()
    assert "Ищем CAS в базе шеринга..." in text

    text = user.get_message_text()
    assert f"Реагентом могут поделиться эти контакты: {user.user.username}" in text

    # у двоих есть
    user.send_message("120-46-7")

    text = user.get_message_text()
    assert "Ищем CAS в базе шеринга..." in text

    text = user.get_message_text()
    assert f"Реагентом могут поделиться эти контакты: {user.user.username}, {admin.user.username}" \
           in text

    user.send_message(CANCEL_SEARCH)
    text = user.get_message_text()
    assert "Поиск завершен" in text


def test_search_smiles(purge_users_collection: None,  # очищаем БД
                       bot, user, admin,
                       mock_pubchempy_get_compounds):

    mock_pubchempy_get_compounds.side_effect = [None, None]

    dbuser = create_user_with_reagents(user, reagents1)
    dbadmin = create_user_with_reagents(admin, reagents2)

    user.send_command("/search")

    text = user.get_message_text()
    assert "Пришли интересующий CAS-номер" in text

    user.send_message("C[C@H]([NH3+])CO")

    text = user.get_message_text()
    assert "Не похоже на CAS. Сейчас поищем по названию..." in text

    text = user.get_message_text()
    assert f"Реагентом могут поделиться эти контакты: {user.user.username}" in text

    # у двоих есть
    user.send_message("O=C(CC(=O)c1ccccc1)c2ccccc2")

    text = user.get_message_text()
    assert "Не похоже на CAS. Сейчас поищем по названию..." in text

    text = user.get_message_text()
    assert f"Реагентом могут поделиться эти контакты: {user.user.username}, {admin.user.username}" \
           in text

    user.send_message(CANCEL_SEARCH)
    text = user.get_message_text()
    assert "Поиск завершен" in text


def test_search_name(purge_users_collection: None,  # очищаем БД
                     bot, user, admin,
                     mock_pubchempy_get_compounds,
                     mock_cirpy_resolve):

    mock_pubchempy_get_compounds.side_effect = [[PubChempyComponent('C1=CC=C(C=C1)C(=O)CC(=O)C2=CC=CC=C2')]]
    mock_cirpy_resolve.side_effect = [['61346-73-4', '120-46-7']]

    dbuser = create_user_with_reagents(user, reagents1)
    dbadmin = create_user_with_reagents(admin, reagents2)

    user.send_command("/search")

    text = user.get_message_text()
    assert "Пришли интересующий CAS-номер" in text

    user.send_message("Dibenzoylmethane")

    text = user.get_message_text()
    assert "Не похоже на CAS. Сейчас поищем по названию..." in text

    text = user.get_message_text()
    assert "Ищем по пользователям SMILES" in text

    text = user.get_message_text()
    assert "Ищем по пользователям CAS=61346-73-4" in text

    text = user.get_message_text()
    assert "Ищем по пользователям CAS=120-46-7" in text

    text = user.get_message_text()
    assert f"Реагентом могут поделиться эти контакты: {user.user.username}, {admin.user.username}" in text

    user.send_message(CANCEL_SEARCH)
    text = user.get_message_text()
    assert "Поиск завершен" in text

