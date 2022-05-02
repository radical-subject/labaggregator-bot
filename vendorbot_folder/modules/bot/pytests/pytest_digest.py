
from io import TextIOWrapper
from .pytest_search import create_user_with_reagents
from modules.db.users import users_collection


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


def get_digests_content(admin, wait_cas_count):
    admin.send_command('/digest')

    text = admin.get_message_text()
    assert 'Ожидайте: список обрабатывается...' in text

    text = admin.get_message_text()
    assert f'Всего {wait_cas_count} CAS' in text

    document = admin.get_document()

    file_io = document.path_or_bytes
    assert file_io.name == 'digest_cas.txt'

    digest_cas = TextIOWrapper(file_io, encoding='utf-8').read()

    document = admin.get_document()

    file_io = document.path_or_bytes
    assert file_io.name == 'digest.txt'

    digest = TextIOWrapper(file_io, encoding='utf-8').read()

    return digest, digest_cas


def test_digest(purge_users_collection: None,  # очищаем БД
                bot, user, admin):

    user.send_command('/digest')

    user.assert_message("digest только для администраторов")

    admin.send_command('/digest')

    text = admin.get_message_text()
    assert 'Ожидайте: список обрабатывается...' in text

    text = admin.get_message_text()
    assert 'Всего 0 CAS' in text

    # теперь будет 3

    dbuser = create_user_with_reagents(user, reagents1)

    digest, digest_cas = get_digests_content(admin, 3)

    assert digest_cas == """120-46-7
2749-11-3
75-64-9"""

    assert digest == """120-46-7 : user1
2749-11-3 : user1
75-64-9 : user1"""

    # будет попрежнему 3, т.к. мы только уникальные выводим
    dbadmin = create_user_with_reagents(admin, reagents1)

    digest, digest_cas = get_digests_content(admin, 3)

    assert digest_cas == """120-46-7
2749-11-3
75-64-9"""

    assert digest == """120-46-7 : user1, user336091411
2749-11-3 : user1, user336091411
75-64-9 : user1, user336091411"""

    # Добавим Админу новые реактивы (новых там 2 шт) итого 5 будет
    dbadmin["user_reagents"].extend(reagents2)
    users_collection.update_user(admin.id, dbadmin)

    digest, digest_cas = get_digests_content(admin, 5)

    assert digest_cas == """106-97-8
120-46-7
13896-65-6
2749-11-3
75-64-9"""

    assert digest == """106-97-8 : user336091411
120-46-7 : user1, user336091411
13896-65-6 : user336091411
2749-11-3 : user1, user336091411
75-64-9 : user1, user336091411"""
