
def test_start_message(bot, user, admin):
    """
    Проверим какие команды показываются обычному юзеру, а какие админу
    """
    user_commands = ['/start', '/help', '/manage', '/search']
    admin_commands = ['/digest', '/dump', '/purge_handler', '/blacklist_reload', '/logs']

    user.send_command('/start')

    message = user.get_message()

    assert message
    answer = message['text']

    assert 'Рады тебя видеть' in answer, answer

    for c in user_commands:
        assert c in answer

    for a in admin_commands:
        assert a not in answer

    #
    admin.send_command('/start')
    message = admin.get_message()

    assert message
    answer = message['text']

    assert 'Рады тебя видеть' in answer, answer

    for c in user_commands:
        assert c in answer

    for a in admin_commands:
        assert a in answer


def test_start_no_username(purge_users_collection: None,  # очищаем БД
                           bot, user, anonim):

    anonim.send_command('/start')

    message = anonim.get_message()
    assert 'Привет' in message['text']

    anonim.assert_get_keyboard('Please share your contact',
                               'Share contact',
                               request_contact=True)

    anonim.send_contact()

    message = anonim.get_message()
    assert 'Thanks for sharing your contact' in message['text']


def test_help(bot, user):

    user.send_command('/help')

    message = user.get_message()

    assert message
    assert 'Добро пожаловать' in message['text']
