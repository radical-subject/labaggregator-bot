

def test_digest(purge_users_collection: None,  # очищаем БД
                bot, user, admin):

    user.send_command('/digest')

    user.assert_message("digest только для администраторов")

    admin.send_command('/digest')

    message = admin.get_message()

    assert message
    assert 'Ожидайте: список обрабатывается...' in message['text']

    message = admin.get_message()

    assert message
    assert 'Всего 0 CAS' in message['text']
