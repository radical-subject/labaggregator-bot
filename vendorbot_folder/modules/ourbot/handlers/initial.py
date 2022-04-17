
from modules.ourbot.logger import logger
from telegram import (ReplyKeyboardMarkup, KeyboardButton, ParseMode)

from telegram import Update
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler, MessageHandler, Filters

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.helpers import bot_commands_text, CONV_START, REQ_CONTACT_STATE
from modules.db.dbmodel import users_collection


class Inital(Handlers):

    """
    класс содержащий в себе стартовые функции хендлеры. наследует класс Handlers,
    в котором прописаны флаги СОСТОЯНИЯ (для диалогов?) и распаковка словаря **db_clients 
    содержащего в себе список подключений к базе данных.
    """

    def __init__(self, bot, db_instances):
        super().__init__(db_instances)
        self.bot = bot

    def start(self, update: Update, context: CallbackContext):
        """
        Стартовая точка общения с ботом.
        welcome message and initialization of user by inserting his data into DB
        """
        user_info = update.message.from_user
        chat_id = update.message.chat.id
        logger.info(f'start({chat_id})')

        # приветственное сообщение юзеру
        text = f"""Привет, {user_info.first_name}! 👩🏻‍💻 
Рады тебя видеть, мхехе.
{bot_commands_text(chat_id)}"""

        update.message.reply_text(text, parse_mode=ParseMode.HTML)

        #
        if not user_info.username:
            logger.info(f'no username({chat_id})')
            reply_markup = ReplyKeyboardMarkup([[KeyboardButton('Share contact', request_contact=True)]])
            self.bot.sendMessage(chat_id, 'You havent set up your username. You will not be able to use sharing. Please share your contact to proceed:', reply_markup=reply_markup)

            return REQ_CONTACT_STATE

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

        return ConversationHandler.END

    def get_contact(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        logger.info(f'get_contact({chat_id})')

        phone_number = update.message.contact.phone_number
        logger.info(phone_number)
        
        reply_markup = ReplyKeyboardMarkup([])
        self.bot.sendMessage(chat_id, 'Thanks.', reply_markup=reply_markup)

        return ConversationHandler.END

    def exit(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        logger.info(f'start.exit({chat_id})')
        return ConversationHandler.END

    def help_command(self, update: Update, context: CallbackContext):
        """Send a message when the command /help is issued."""
        chat_id = update.message.chat_id
        logger.info(f'help({chat_id})')

        update.message.reply_text("""
Добро пожаловать в альфа-версию бота для обмена реактивов. 
Чтобы получить доступ к системе обмена необходимо поделиться своим списком. 
/manage - Загрузить свой список реагентов можно в виде .txt файла с CAS номерами в столбик. 
/search - Поиск по CAS по базе реагентов присланных для обмена.

Общественные списки также публикуются дайджестами в канале Лабаггрегатора @labaggregator.
        """, parse_mode=ParseMode.HTML)

    def register_handler(self, dispatcher):
        
        dispatcher.add_handler(CommandHandler('help', self.help_command))

        self.conversation_handler = ConversationHandler(
            entry_points=[CommandHandler('start', self.start)],
            states={
                REQ_CONTACT_STATE: [
                    MessageHandler(Filters.contact, self.get_contact)
                ],
            },
            fallbacks=[MessageHandler(Filters.command, self.exit),
                       MessageHandler(Filters.text, self.exit)],
        )
        
        dispatcher.add_handler(self.conversation_handler, CONV_START)
