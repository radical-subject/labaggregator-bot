import io
from typing import List, Optional, Tuple, Any

import pandas as pd
from modules.reagent import Reagent
from telegram import Update
from telegram.ext import CallbackContext

# у каждого conversational handler должна быть своя группа
CONV_START, CONV_MANAGE, CONV_SEARCH, CONV_APPEND = range(4)

# TODO у всех conversational handler состояний должен быть свой номер.
# кажется проблема с асинхронными хэндлерами
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
/logs - скачать bot.log
/append - добавить новые компоненты в список
/digest - загрузить все shared
/purge_handler - очистка бд (только админам)
/dump - дамп базы данных (присылает в лс зип-дамп)
/blacklist_reload - заполнение базы блеклиста и обсчет. команда выполняется асинхронно
"""
    return text


LIST_OF_ADMINS = [
    336091411,
    122267418,
    588927967,
    47390523,
    250291302,
    880726373,
    490796163,
]  # tg id's 336091411, 122267418, 588927967, 47390523


def is_admin_chat(chat_id):
    return chat_id in LIST_OF_ADMINS


def get_file(
    update: Update, context: CallbackContext
) -> [Optional[io.BytesIO], Optional[str]]:
    """
    Читаем файл, отправленный в бот
    :param update:
    :param context:
    :return:
    """
    document = update.message.document
    if document:
        out = io.BytesIO()
        context.bot.get_file(document).download(out=out)
        out.seek(0)
        return out, document["file_name"]
    return None, None


def file_to_dataframe(file: io.BytesIO, name: str) -> Optional[pd.DataFrame]:
    """
    Читаем текстовый файл пандой
    :param file: файл
    :param name: имя файла
    :return: pd.DataFrame
    """
    ext = name.split(".")[-1]
    if ext == "txt":
        return pd.read_csv(file, names=["CAS"], usecols=[0])
    elif ext == "csv":
        return pd.read_csv(
            file,
            delimiter=";|,|\t",
            usecols=["CAS", "location", "name"],
            encoding="unicode_escape",
        )
    elif ext in ["xlsx", "xls"]:
        return pd.read_excel(file, usecols=["CAS", "location", "name"])


def get_contact_from_dataframe(df: pd.DataFrame) -> Tuple[Optional[str], pd.DataFrame]:
    """
    # оставляю возможность хардкодить вручную контакт,
    прописывая первую строку импортируемого файла руками:
    # в формате reagents_contact:+79265776746
    :param df: содержимое файла
    :param user_info:
    :return:
    """
    try:
        if "CAS" in df:
            contact = df["CAS"][0]
            if contact.startswith("reagents_contact:"):
                df = df.iloc[
                    1:
                ]  # удаляем 1ю строчку  TODO а надо ли ? мы уже считали же всё
            return contact.split(":")[1], df
    except:
        return None, df


def parse_cas(df_value: Any) -> str:
    """
    При пустой ячейке там может прочитаться float!
    :param df_value: ячейка excel
    :return: строка
    """
    if isinstance(df_value, str):
        return df_value.split("\n")[0]


def df_to_reagents(df: pd.DataFrame) -> List[Reagent]:
    """
    TODO: мы получили из эксель 100 строк
    но только 90 строк признали хорошими
    юзеру будет интересно какие строки признаны плохими,
    кажется стоит вернуть ему эксель файл или написать номер строк что ли..
    :param df: прочитанный файл
    :return: список Reagent
    """
    out = []
    for index, row in df.iterrows():
        r = Reagent()
        if 'CAS' in row.index:
            r.cas = parse_cas(row['CAS'])
            if not r.cas:
                continue   # если нет CAS, то не вставляем в БД

        if 'location' in row.index:
            r.location = str(row['location'])
        if 'name' in row.index:
            r.name = str(row['name'])

        out.append(r)
    return out
