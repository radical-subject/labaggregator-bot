

def test_digest(purge_users_collection: None,  # очищаем БД
                bot, user, admin):
    user.init_dialog()
    admin.init_dialog()

    user.send_command('/digest')

    message = user.get_message()
    assert not message, f'/digest только для администраторов, {message}'

    admin.send_command('/digest')

    message = admin.get_message()

    assert message
    assert 'Ожидайте: список обрабатывается...' in message['text']

    message = admin.get_message()

    assert message
    assert 'Всего 0 CAS' in message['text']
