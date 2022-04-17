
from telegram import Update, ReplyKeyboardMarkup, KeyboardButton, ReplyKeyboardRemove, \
    InlineKeyboardMarkup, InlineKeyboardButton, ParseMode
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler, \
    RegexHandler, MessageHandler, CallbackQueryHandler, Filters

from modules.ourbot.handlers.helpers import CONV_SEARCH
from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.helpers import is_CAS_number
from modules.ourbot.logger import logger

from modules.db import dbschema
from modules.db.dbmodel import users_collection


SEARCH_STATE = range(1)
CANCEL_CALLBACK = str('SEARCH:CANCEL')

cancel_keyboard = [
    [
        InlineKeyboardButton("CANCEL SEARCH", callback_data=CANCEL_CALLBACK)
    ]
]


class Search(Handlers):
    def __init__(self, bot, db_instances):
        super(Search, self).__init__(db_instances)
        self.collection = "users_collection"

    def search(self, update: Update, context: CallbackContext):
        """
        Старт ветки диалога "поиск"
        """
        chat_id = update.message.chat_id
        logger.info(f'search({chat_id})')
        reply_markup = InlineKeyboardMarkup(cancel_keyboard)
        update.message.reply_text("🙋🏻‍♀️ Enter query (name or CAS):\n\n"
                                  "🖋 Пришли интересующий CAS-номер:",
                                  reply_markup=reply_markup)
        return SEARCH_STATE

    def search_cas(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        logger.info(f'search_cas({chat_id})')

        text = update.message.text
        if is_CAS_number(text):
            update.message.reply_text('Ищем CAS в базе шеринга...')

            users = users_collection.get_users_by_cas(text)

            contacts = []
            for user in users:
                user_reagents_object = dbschema.UserReagents(**user)

                for contact in user_reagents_object.get_contacts_for_CAS(text):
                    if contact not in contacts:
                        contacts.append(contact)

            if contacts:
                update.message.reply_text(f'Реагентом могут поделиться эти контакты: {", ".join(contacts)}')
            else:
                update.message.reply_text('Реагентом пока никто не готов поделиться.')

        else:
            update.message.reply_text('Неправильный CAS номер. Попробуйте еще раз.')

        return SEARCH_STATE

    def exit(self, update: Update, context: CallbackContext) -> int:
        """
        Выход из ветки диалога "поиск"
        """
        chat_id = update.message.chat_id
        logger.info(f'search.exit({chat_id})')

        # необходимо согласно мануалу ответить на query
        query = update.callback_query
        query.answer()

        # берем последнее сообщение бота
        sent_message = update.callback_query.message

        # редактируем его меняя текст и убирая кнопку. диалог завершен.
        context.bot.edit_message_text(
            text=f'STOPPED',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=None,
            parse_mode=ParseMode.MARKDOWN
        )
        
        # now clear all cached data
        # clear assosiated with user data and custom context variables
        context.chat_data.clear()
        context.user_data.clear()

        return ConversationHandler.END

    def register_handler(self, dispatcher):

        conv_search = ConversationHandler(
            entry_points=[CommandHandler('search', self.search),],
            states={
                SEARCH_STATE: [
                    CallbackQueryHandler(self.exit, pattern=CANCEL_CALLBACK),
                    MessageHandler(Filters.text & ~Filters.command, self.search_cas, run_async=True)
                ],
            },
            fallbacks=[MessageHandler(Filters.command, self.exit),
                       MessageHandler(Filters.text, self.exit)],
        )

        dispatcher.add_handler(conv_search, CONV_SEARCH)
