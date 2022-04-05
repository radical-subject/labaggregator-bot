
from telegram import Update, ReplyKeyboardMarkup, KeyboardButton, ReplyKeyboardRemove
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler, \
    RegexHandler, MessageHandler, Filters

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.helpers import is_CAS_number
from modules.ourbot.logger import log

from modules.db import dbmodel, dbschema

SEARCH_STATE = range(1)
CANCEL = 'Отмена'
CANCEL_REGEXP = '^Отмена$'

cancel_keyboard = [[KeyboardButton(CANCEL)]]


class Search(Handlers):
    def __init__(self, bot, db_instances):
        super(Search, self).__init__(db_instances)
        self.bot = bot
        self.collection = "users_collection"

    def search(self, update: Update, context: CallbackContext):
        """
        Старт ветки диалога "поиск"
        """
        update.message.reply_text("🙋🏻‍♀️ Enter query (name or CAS):\n\n"
                                  "🖋 Пришли интересующий CAS-номер:",
                                  reply_markup=ReplyKeyboardMarkup(cancel_keyboard, resize_keyboard=True))
        return SEARCH_STATE

    def search_cas(self, update: Update, context: CallbackContext):

        chat_id = update.message.chat_id

        text = update.message.text
        if is_CAS_number(text):
            update.message.reply_text('ищем CAS в базе шеринга...')
            mongo_query = {"user_reagents": { '$elemMatch': { 'CAS': text}}}
            result = dbmodel.get_records(self.vendorbot_db_client, self.db_instances["vendorbot_db"], self.collection, mongo_query)
            try:
                update.message.reply_text(f'Реагентом могут поделиться эти контакты: {result[0]["username"]}')
            except AttributeError: 
                update.message.reply_text('Реагентом пока никто не готов поделиться.')

        else:
            update.message.reply_text('Неправильный CAS номер. Попробуйте еще раз.')

        return SEARCH_STATE

    def cancel(self, update: Update, context: CallbackContext) -> int:
        """
        Выход из ветки диалога "поиск"
        """
        log.info(f"Завершение поиска")
        update.message.reply_text("Завершение поиска",
                                  reply_markup=ReplyKeyboardRemove())

        return ConversationHandler.END

    def register_handler(self, dispatcher):

        self.conversation_handler = ConversationHandler(
            entry_points=[CommandHandler('search', self.search),],
            states={
                SEARCH_STATE: [MessageHandler(Filters.regex(CANCEL_REGEXP), self.cancel),
                               MessageHandler(Filters.text & ~Filters.command, self.search_cas,
                                              run_async=True)
                               ],
            },
            fallbacks=[RegexHandler(CANCEL_REGEXP, self.cancel)],
        )

        dispatcher.add_handler(self.conversation_handler, 1)
