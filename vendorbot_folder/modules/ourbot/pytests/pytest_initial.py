
def test_start(bot, user, admin):
    """
    Проверим какие команды показываются обычному юзеру, а какие админу
    """

    user_commands = ['/start', '/help', '/manage', '/search']
    admin_commands = ['/digest', '/dump', '/purge_handler', '/blacklist_update']

    user.send_command('/start')

    message = user.get_message()

    assert message
    answer = message['text']

    assert 'Рады тебя видеть, мхехе.' in answer, answer

    for c in user_commands:
        assert c in answer

    for a in admin_commands:
        assert a not in answer

    #
    admin.send_command('/start')

    message = admin.get_message()

    assert message
    answer = message['text']

    assert 'Рады тебя видеть, мхехе.' in answer, answer

    for c in user_commands:
        assert c in answer

    for a in admin_commands:
        assert a in answer


def test_help(bot, user):

    user.send_command('/help')

    message = user.get_message()

    assert message
    assert 'Добро пожаловать' in message['text']
