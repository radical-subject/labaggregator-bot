
from typing import List

import io
from csv import reader, excel
from telegram.ext import CallbackContext
from telegram import Update

# у каждого conversational handler должна быть своя группа
CONV_START, CONV_MANAGE, CONV_SEARCH, CONV_APPEND = range(4)

# у всех conversational handler состояний должен быть свой номер. кажется проблема с асинхронными хэндлерами
UPLOAD_STATE, REQ_CONTACT_STATE, SEARCH_STATE, APPEND_STATE = range(4)


def bot_commands_text(chat_id):
    text = """
Доступны следующие команды:
/start - приветствие
/help - инструкции по пользованию
/manage - загрузить список компонентов (заменит старый)
/search - поиск по CAS
"""
    if is_admin_chat(chat_id):
        text += """\n== Админам ==
/append - добавить новые компоненты в список
/digest - загрузить все shared
/purge_handler - очистка бд (только админам)
/dump - дамп базы данных (присылает в лс зип-дамп)
/blacklist_update - заполнение базы блеклиста и обсчет. команда выполняется асинхронно
"""
    return text


LIST_OF_ADMINS = [336091411, 122267418, 588927967, 47390523, 250291302]  # tg id's 336091411, 122267418, 588927967, 47390523


def is_admin_chat(chat_id):
    return chat_id in LIST_OF_ADMINS


def get_csv_content(update: Update, context: CallbackContext) -> List[List[str]]:
    """
    Возвращает содержимое csv файла присланного пользователем
    :param update:
    :param context:
    :return:
    """
    out = io.BytesIO()
    context.bot.get_file(update.message.document).download(out=out)
    out.seek(0)

    content = out.read().decode('utf-8')
    rows = []
    c = reader(content.splitlines(), delimiter=',')  # out.read().decode('utf-8').splitlines()
    for r in c:
        rows.append(r)
    return rows


def get_txt_content(update: Update, context: CallbackContext) -> List[str]:
    out = io.BytesIO()
    context.bot.get_file(update.message.document).download(out=out)
    out.seek(0)

    content = out.read().decode('utf-8')
    rows = []
    for r in content.splitlines():
        rows.append(r.strip())
    return rows
