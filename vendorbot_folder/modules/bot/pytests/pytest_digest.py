
from .pytest_search import create_user_with_reagents

reagents1 = [{'reagent_internal_id': '1', 'CAS': '2749-11-3', 'SMILES': 'C[C@H]([NH3+])CO',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'},
             {'reagent_internal_id': '2', 'CAS': '75-64-9', 'SMILES': 'CC(C)(C)N',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'},
             {'reagent_internal_id': '3', 'CAS': '120-46-7', 'SMILES': 'O=C(CC(=O)c1ccccc1)c2ccccc2',
              'sharing_status': 'shared', 'timestamp': '22.04.2022 14:26'}]


def test_digest(purge_users_collection: None,  # очищаем БД
                bot, user, admin):

    user.send_command('/digest')

    user.assert_message("digest только для администраторов")

    admin.send_command('/digest')

    text = admin.get_message_text()
    assert 'Ожидайте: список обрабатывается...' in text

    text = admin.get_message_text()
    assert 'Всего 0 CAS' in text

    """
    НУЖНО реализовать отправку файлов ботом
    # теперь будет 3
    dbuser = create_user_with_reagents(user, reagents1)

    admin.send_command('/digest')

    text = admin.get_message_text()
    assert 'Ожидайте: список обрабатывается...' in text

    text = admin.get_message_text()
    assert 'Всего 3 CAS' in text

    # теперь будет 6
    dbadmin = create_user_with_reagents(admin, reagents1)

    admin.send_command('/digest')

    text = admin.get_message_text()
    assert 'Ожидайте: список обрабатывается...' in text

    text = admin.get_message_text()
    assert 'Всего 6 CAS' in text
    """
