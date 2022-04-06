
import json
from bson import ObjectId
import pandas as pd

from modules.ourbot.logger import logger
from telegram import (ReplyKeyboardMarkup, KeyboardButton, ParseMode)

from telegram import Update, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler

from modules.db.dbschema import UserReagents
        
from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.decorators import log_errors
from modules.ourbot.handlers.helpers import is_admin_chat
from modules.db import dbmodel, dbschema


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
        self.bot=bot
        self.collection = "users_collection"

    @log_errors
    def start(self, update: Update, context: CallbackContext):
        """
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

        if not dbmodel.get_user(user_info.id):
            dbmodel.add_user(userdata)
        else:
            logger.info("User already exists: skipping insertion of userdata in DB")

        # associated with user chat and context stored data should be cleaned up to prevent mess
        context.chat_data.clear()
        context.user_data.clear()

        return self.INITIAL  #todo не нужно - это не conversation_handler !

    #TODO remove
    def exit_command(self, update: Update, context: CallbackContext):
        """
        hadler for terminating all dialog sequences
        """
        try:
            query = update.callback_query
            if query != None:
                reply_markup = InlineKeyboardMarkup([])
                query.edit_message_text(
                    text="You cancelled db removal. Да поможет тебе святой Януарий!",
                    reply_markup=reply_markup
                )
            else:
                update.message.reply_text(f"""Выход из диалога. Да поможет тебе святой Антоний.""")
        except:
            update.message.reply_text(f"""Выход из диалога. Да поможет тебе святой Антоний.""")
            pass
        # now clear all cached data
        # clear assosiated with user data and custom context variables
        context.chat_data.clear()
        context.user_data.clear()
        return ConversationHandler.END #todo не нужно - это не conversation_handler !

    #TODO remove
    def my_lab(self, update: Update, context: CallbackContext):
        current_lab = context.user_data.get('current_lab')
        result = 'None' if current_lab is None else JSONEncoder().encode(current_lab)
        update.message.reply_text(result, parse_mode='HTML')

    def help_command(self, update: Update, context: CallbackContext):
        """Send a message when the command /help is issued."""
        update.message.reply_text("""
Добро пожаловать в альфа-версию бота для обмена реактивов. 
Чтобы получить доступ к системе обмена необходимо поделиться своим списком. 
/manage - Загрузить свой список реагентов можно в виде .txt файла с CAS номерами в столбик. 
/search - В разработке. Текстовый поиск произвольного формата по базе общественных реагентов.

пока /search в разработке, общественные списки будут публиковаться дайджестами в канале лабаггрегатора.
        """, parse_mode=ParseMode.HTML)

        return self.INITIAL #todo не нужно - это не conversation_handler !

    @log_errors
    # @run_async
    def resolve_tests(self, update: Update, context: CallbackContext):

        input_txt_file_path = "./srs/user_reagent_lists_import/Chusov_1.txt"
        import_CAS_df = pd.read_csv(input_txt_file_path, header = None)
        CAS_list = import_CAS_df[0].tolist()

        # retrieving data from user message
        # ищем запись относящуюся к пользователю
        user_id = update.message.from_user.id
        mongo_query = {"user_id": user_id}
        logger.info(f"mongo_query: {mongo_query}")
        user_info = update.message.from_user
        chat_id = update.message.chat.id

        initial_record = dbmodel.get_records(self.vendorbot_db_client, self.db_instances["vendorbot_db"], self.collection, mongo_query)
        logger.info(f"initial_record: {initial_record}")
        user_reagents_object = UserReagents(**initial_record[0])
        
        # импорт листа реагентов с фильтрациями
        user_reagents_object.add_list_of_reagents(user_info.id, user_info.username, self.blacklist_rdkit_db_client, self.db_instances["blacklist_rdkit_db"], CAS_list)
        # экспорт JSON - не работает с pymongo! нужен dict
        data = user_reagents_object.export()
        # записываем в базу объект 
        dbmodel.update_record(self.vendorbot_db_client, self.db_instances["vendorbot_db"], self.collection, mongo_query, data)
        update.message.reply_text(f"{user_reagents_object.shared_reagents()[0]}")

    @log_errors
    def capture_contact(self, update: Update, context: CallbackContext):
        # retrieving data from user message
        user_info = update.message.from_user
        chat_id = update.message.chat.id

        reply_markup = ReplyKeyboardMarkup([[KeyboardButton('Share contact', request_contact=True)]])
        self.bot.sendMessage(chat_id, 'Example', reply_markup=reply_markup)

    @log_errors
    def set_tag(self, update: Update, context: CallbackContext):
        """set tag for timerdata"""
        pass

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('start', self.start))
        dispatcher.add_handler(CommandHandler('my_lab', self.my_lab))
        dispatcher.add_handler(CommandHandler('end', self.exit_command))
        dispatcher.add_handler(CommandHandler('help', self.help_command))
        dispatcher.add_handler(CommandHandler('resolve_tests', self.resolve_tests))
        dispatcher.add_handler(CommandHandler('capture_contact', self.capture_contact))
