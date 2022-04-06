
from typing import List

import io
from csv import reader, excel
from telegram.ext import CallbackContext
from telegram import Update


#TODO don't use
def remove_job_if_exists(name: str, context: CallbackContext) -> bool:
    """
    Удаляет задачу из очереди задач Бота по имени

    :param name: имя задачи. должно быть уникальным.
    :param context:
    :return: False - нет такой задачи, True - задача удалена
    """
    current_jobs = context.job_queue.get_jobs_by_name(name)
    if not current_jobs:
        return False
    for job in current_jobs:
        job.schedule_removal()
    return True


LIST_OF_ADMINS = [336091411, 122267418, 588927967, 47390523]  # tg id's 336091411, 122267418, 588927967, 47390523


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
