
import logging
logger = logging.getLogger(__name__)
from telegram import (ReplyKeyboardMarkup, KeyboardButton, ParseMode)

from telegram import Update
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler, MessageHandler, Filters

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.helpers import bot_commands_text, CONV_START, REQ_CONTACT_STATE
from modules.db.dbmodel import users_collection
from modules.db.dbschema import UserReagents


def user_from_user_info(user_info, phone_number: str = ""):
    return {
        "_id": user_info.id,
        "user_id": user_info.id,
        "username": user_info.username,
        "firstname": user_info.first_name,
        "lastname": user_info.last_name,
        "phone_number": phone_number
    }


class Initial(Handlers):

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
Рады тебя видеть. Этот бот помогает ученым делиться друг с другом образцами химреактивов.
{bot_commands_text(chat_id)}"""

        update.message.reply_text(text, parse_mode=ParseMode.HTML)

        user = users_collection.get_user(user_info.id)
        if not user:
            # запись данных юзера в БД произойдет сразу, если у юзера есть юзернейм
            user = user_from_user_info(user_info)
            users_collection.add_user(user)
            logger.info(f"user added")
        else:
            logger.info(f"User ({user_info.id}) exists")

        logger.info(f"userdata({user})")

        # проверка на наличие юзернейма, если его не предоставлено - идет запрос контакта (телефонного номера)
        if not user["username"] and ("phone_number" not in user or not user["phone_number"]):
            logger.info(f"no username or no phone ({chat_id})")

            reply_markup = ReplyKeyboardMarkup([[KeyboardButton("Share contact", request_contact=True)]],
                                               resize_keyboard=True, one_time_keyboard=True)

            self.bot.sendMessage(chat_id, "You haven\'t setup your username. You will not be able to use sharing. "
                                          "Please share your contact to proceed any further:",
                                 reply_markup=reply_markup)
            return REQ_CONTACT_STATE

        # associated with user chat and context stored data should be cleaned up to prevent mess
        context.chat_data.clear()
        context.user_data.clear()

        return ConversationHandler.END

    def get_contact(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        user_info = update.message.from_user
        logger.info(f"get_contact({chat_id})")

        phone_number = update.message.contact.phone_number
        logger.info(f"get_contact phone={phone_number}")

        user = users_collection.get_user(user_info.id)
        if not user:
            logger.error("user должен был быть создан в /start. проверить!")
            user = user_from_user_info(user_info, phone_number=phone_number)
            users_collection.add_user(user)
        else:
            user["phone_number"] = phone_number
            users_collection.update_user(user_info.id, user)

        self.bot.sendMessage(chat_id, "Thanks for sharing your contact. "
                                      "Now you will be able to upload your list of reagents. /manage")

        context.chat_data.clear()
        context.user_data.clear()
        return ConversationHandler.END

    def exit(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        logger.info(f"start.exit({chat_id})")

        context.chat_data.clear()
        context.user_data.clear()
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
