
from typing import List

import io
#from csv import reader
from telegram.ext import CallbackContext
from telegram import Update

from pandas import read_excel, read_csv, DataFrame

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
/calculate_hashes - подготовить базу к поиску
/append - добавить новые компоненты в список
/digest - загрузить все shared
/purge_handler - очистка бд (только админам)
/dump - дамп базы данных (присылает в лс зип-дамп)
/blacklist_reload - заполнение базы блеклиста и обсчет. команда выполняется асинхронно
"""
    return text


LIST_OF_ADMINS = [336091411, 122267418, 588927967, 47390523, 250291302, 880726373, 490796163]  # tg id's 336091411, 122267418, 588927967, 47390523


def is_admin_chat(chat_id):
    return chat_id in LIST_OF_ADMINS


def get_file_content(update: Update, context: CallbackContext) -> List[List[str]]:
    out = io.BytesIO()
    document = update.message.document
    context.bot.get_file(document).download(out=out)
    out.seek(0)

    ext = document['file_name'].split('.')[-1]
    if ext == 'txt':
        return read_csv(out, names=['CAS'], usecols=[0])
    elif ext == 'csv':
        return read_csv(out, delimiter=';|,|\t', usecols=['CAS', 'location', 'name'], encoding= 'unicode_escape')
    elif ext == 'xlsx':
        return read_excel(out, usecols=['CAS', 'location', 'name'])
    else:
        return 'Valid file extensions are [ .txt | .csv | .xlsx ]'