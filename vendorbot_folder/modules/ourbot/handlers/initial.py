import logging
import pymongo
from telegram import Update, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors
from modules.db import dbmodel, dbschema
import json
from bson import ObjectId

logger = logging.getLogger(__name__)

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

    def __init__(self, db_instances):
        super().__init__(db_instances)
        self.collection = "users_collection"
        self.collection_2 = "timer_data_collection"

    @log_errors
    def start_msg(self, update: Update, context: CallbackContext):
        """
        welcome message and initialization of user by inserting his data into DB
        """
        # retrieving data from user message
        user_info = update.message.from_user
        chat_id = update.message.chat.id

        # приветственное сообщение юзеру
        update.message.reply_text(
            """Привет, {}! 👩🏻‍💻
Доступны следующие команды:
/start - приветствие 
/purge_handler - очистка бд (только админам)
/help - инструкции по пользованию
/dump - дамп базы данных (присылает в лс зип-дамп)
/timer - запуск диалога таймера
/choose_category - то же самое что выбор лаборатории, только пока не доделано до конца.
/edit_category - архивирование, разархивирование, удаление, (добавить переименование)
/today - общая информация за сегодня о проделанных делах, время нетто, время брутто (заглушка == 10 часов в день), список записей таймера
==================================
/choose_lab - выбираем лабораторию
/my_lab - показывает выбранную лабораторию. осторожно с переменными состояния хранящимися в контексте - там уже адская путаница. к тому же некоторые команды норовят переменные состояния сбросить .clear(). """.format(user_info.first_name), parse_mode='HTML')

        # запись данных юзера в БД
        userdata_dict = {
            "_id": user_info.id,
            "user_id": user_info.id,
            "username": "@{}".format(user_info.username),
            "firstname": user_info.first_name,
            "lastname": user_info.last_name
        }
        try:
            # logger.info(f"{self.timerbot_db_client}, {self.db_instances['timerbot_db']}, {self.collection}, {userdata_dict}")
            dbmodel.add_records(self.root_client, self.db_instances["timerbot_db"], self.collection, userdata_dict)
            logger.info('user initialized by /start command.')
        except pymongo.errors.DuplicateKeyError:
            logger.info("User already exists: skipping insertion of userdata in DB")
        # associated with user chat and context stored data should be cleaned up to prevent mess
        context.chat_data.clear()
        user_data = context.user_data
        user_data.clear()

        return self.INITIAL

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
        # equivalent of return ConversationHandler.END
        return -1

    def my_lab(self, update: Update, context: CallbackContext):
        current_lab = context.user_data.get('current_lab')
        result = 'None' if current_lab is None else JSONEncoder().encode(current_lab)
        update.message.reply_text(result, parse_mode='HTML')

    def help_command(self, update: Update, context: CallbackContext):
        """Send a message when the command /help is issued."""
        update.message.reply_text(
        """
При начале диалога таймера таймер запускается сразу же.
В зависимости от того, выбрана ли категория изначально, в конце работы таймера 
формируется объект с информацией о:
1) израсходованном времени (обязательно), 
2) категории (опционально), 
3) названии конкретного дела (обязательно) - заполняется после остановки таймера, 
4) результатах (комментарий, опционально) - заполняется после остановки таймера,
5) timestamp (обязательно). 

любая другая команда абортирует диалог таймера, а любое текстовое сообщение во время работы таймера фильтруется и игнорируется ботом.
        """,
        parse_mode='HTML'
        )
        return self.INITIAL


    @log_errors
    def today_stats(self, update: Update, context: CallbackContext):
        """Send information about work entries for today"""
        user_info = update.message.from_user
        user_id = user_info.id
        # ищем запись относящуюся к пользователю
        mongo_query = {"user_id": user_id}
        previous_records=dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection_2, mongo_query)
        # logger.info(previous_records[0])
        timer_object = dbschema.TimerData(
            **previous_records[0]
        )
        # logger.info(timer_object.export())

        try:
            data = timer_object.export()
        except:
            data = "No records"
        try:
            time_netto_today = timer_object.get_netto_today()
            logger.info(f"time netto today = {timer_object.get_netto_today()} minutes")
        except:
            time_netto_today = 0

        # пересчет в красивый формат
        time_netto_today_hours = time_netto_today // 60
        time_netto_today_minutes = time_netto_today % 60
        logger.info(data)

        str_line = ''
        for i in (list(item.items()) for item in data['timerdata']):
            str_line += "\n".join([': '.join(map(str, tuple_item)) for tuple_item in i])
            str_line += "\n\n" 
        entries = str_line.rstrip("\n")

        update.message.reply_text(f"===================\n{entries}\n===================\n\nTime brutto today == 10 hours.\nthis is temporarily hardcoded.\n\ntime_netto_today = {time_netto_today_hours:.0f} h. {time_netto_today_minutes:.2f} min.")
        
        return self.INITIAL


    @log_errors
    def set_tag(self, update: Update, context: CallbackContext):
        """set tag for timerdata"""
        pass

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('start', self.start_msg))
        dispatcher.add_handler(CommandHandler('my_lab', self.my_lab))
        dispatcher.add_handler(CommandHandler('end', self.exit_command))
        dispatcher.add_handler(CommandHandler('help', self.help_command))
        dispatcher.add_handler(CommandHandler('today', self.today_stats))
