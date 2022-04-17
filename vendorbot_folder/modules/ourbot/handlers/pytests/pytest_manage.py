import os

FILES_DIR = os.path.join(os.path.dirname(__file__), 'files')


def test_manage(purge_users_collection: None,  # очищаем БД
                bot, user, admin):

    user.init_dialog()
    admin.init_dialog()
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

