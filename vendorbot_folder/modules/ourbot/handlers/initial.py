
import json
from bson import ObjectId
import pandas as pd

from modules.ourbot.logger import logger
from telegram import (ReplyKeyboardMarkup, KeyboardButton, ParseMode)

from telegram import Update, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler

from modules.db.dbschema import UserReagents
        
from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.helpers import is_admin_chat
from modules.db.dbmodel import users_collection


class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)
        return json.JSONEncoder.default(self, o)


class Inital(Handlers):

    """
    класс содержащий в себе стартовые функции хендлеры. наследует класс Handlers,
    в котором прописаны флаги СОСТОЯНИЯ (для диалогов?) и распаковка словаря **db_clients 
    содержащего в себе список подключений к базе данных.
    """

    def __init__(self, bot, db_instances):
        super().__init__(db_instances)

    def start(self, update: Update, context: CallbackContext):
        """
        Стартовая точка общения с ботом.
        welcome message and initialization of user by inserting his data into DB
        """
        user_info = update.message.from_user
        chat_id = update.message.chat.id

        # приветственное сообщение юзеру
        text = f"""Привет, {user_info.first_name}! 👩🏻‍💻 
Рады тебя видеть, мхехе.
Доступны следующие команды:
/start - приветствие
/help - инструкции по пользованию
/manage - загрузить список компонентов для обмена
/search - поиск 
"""
        if is_admin_chat(chat_id):
            text += """\n== Админам ==
/digest - загрузить все shared
/purge_handler - очистка бд (только админам)
/dump - дамп базы данных (присылает в лс зип-дамп)
/blacklist_update - заполнение базы блеклиста и обсчет. команда выполняется асинхронно
"""
        update.message.reply_text(text, parse_mode=ParseMode.HTML)

        #TODO А давай писать, если не заполнено user_info.username, то и поиском пользоваться нельзя?
        # иначе мы только по номеру диалога будем знать чей реактив, либо по номеру мобилки

        # запись данных юзера в БД
        userdata = {
            "_id": user_info.id,
            "user_id": user_info.id,
            "username": "@{}".format(user_info.username),
            "firstname": user_info.first_name,
            "lastname": user_info.last_name
        }

        if not users_collection.get_user(user_info.id):
            users_collection.add_user(userdata)
        else:
            logger.info("User already exists: skipping insertion of userdata in DB")

        # associated with user chat and context stored data should be cleaned up to prevent mess
        context.chat_data.clear()
        context.user_data.clear()

    def help_command(self, update: Update, context: CallbackContext):
        """Send a message when the command /help is issued."""
        update.message.reply_text("""
Добро пожаловать в альфа-версию бота для обмена реактивов. 
Чтобы получить доступ к системе обмена необходимо поделиться своим списком. 
/manage - Загрузить свой список реагентов можно в виде .txt файла с CAS номерами в столбик. 
/search - В разработке. Текстовый поиск произвольного формата по базе общественных реагентов.

пока /search в разработке, общественные списки будут публиковаться дайджестами в канале лабаггрегатора.
        """, parse_mode=ParseMode.HTML)

    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('start', self.start))
        dispatcher.add_handler(CommandHandler('help', self.help_command))
