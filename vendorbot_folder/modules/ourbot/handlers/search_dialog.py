import os
import traceback
from telegram import Update, ReplyKeyboardMarkup, KeyboardButton, ReplyKeyboardRemove
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler, \
    MessageHandler, Filters

from modules.ourbot.handlers.helpers import CONV_SEARCH, SEARCH_STATE
from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.cas_to_smiles import what_reagent
import logging

from modules.db.dbschema import get_reagent_contacts
from modules.db.dbmodel import users_collection

logger = logging.getLogger(__name__)

CANCEL_SEARCH = 'Завершить поиск'
cancel_keyboard = [[KeyboardButton(CANCEL_SEARCH)]]


DBSIZE_OPEN_SEARCH = int(os.getenv('DBSIZE_OPEN_SEARCH', 10))


class Search(Handlers):
    def __init__(self, bot, db_instances):
        super(Search, self).__init__(db_instances)
        self.collection = "users_collection"

    def search(self, update: Update, context: CallbackContext):
        """
        Старт ветки диалога "поиск"
        """
        chat_id = update.message.chat_id
        logger.info(f"search({chat_id})")

        user_id = update.message.from_user.id
        count = len(users_collection.get_reagents(user_id))

        if count < DBSIZE_OPEN_SEARCH:  #  and not is_admin_chat(chat_id)
            update.message.reply_text(f"Чтобы разблокировать шеринг, вам необходимо загрузить "
                                      f"в базу не менее {DBSIZE_OPEN_SEARCH} ваших позиций. /manage")
            return ConversationHandler.END

        reply_markup = ReplyKeyboardMarkup(cancel_keyboard, resize_keyboard=True)
        update.message.reply_text("🙋🏻‍♀️ Enter query (name or CAS):\n\n"
                                  "🖋 Пришли интересующий CAS-номер:",
                                  reply_markup=reply_markup)

        return SEARCH_STATE

    def search_cas(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        text = update.message.text

        logger.info(f"search_cas({chat_id}): {text}")

        try:
            contacts = []

            cas_list, smiles_list = what_reagent(text)

            text = "Ищем по пользователям:"
            if cas_list:
                text += f"\nCAS: {' ,'.join(cas_list)}"
            if smiles_list:
                text += f"\nSMILES: {' ,'.join(smiles_list)}"
            update.message.reply_text(text)

            for cas in cas_list:
                contacts.extend(get_reagent_contacts(users_collection.get_users_by_cas(cas), cas))

            for smiles in smiles_list:
                contacts.extend(get_reagent_contacts(users_collection.get_users_by_smiles(smiles), smiles))

            contacts = list(set(contacts))
            if contacts:
                update.message.reply_text(f"Реагентом могут поделиться эти контакты: {', '.join(contacts)}")
            else:
                update.message.reply_text("Реагентом пока никто не готов поделиться.")

        except Exception as err:
            logger.error(traceback.format_exc())
            update.message.reply_text("Ошибка поиска. Похвастайтесь админу, что сломали бот.")

        return SEARCH_STATE

    def exit(self, update: Update, context: CallbackContext) -> int:

        chat_id = update.message.chat_id
        logger.info(f'search.exit({chat_id})')

        update.message.reply_text("Поиск завершен",
                                  reply_markup=ReplyKeyboardRemove())

        return ConversationHandler.END

    def register_handler(self, dispatcher):

        conv_search = ConversationHandler(
            entry_points=[CommandHandler('search', self.search),],
            states={
                SEARCH_STATE: [
                    MessageHandler(Filters.regex(CANCEL_SEARCH), self.exit),
                    MessageHandler(Filters.text & ~Filters.command, self.search_cas, run_async=True)
                ],
            },
            fallbacks=[MessageHandler(Filters.command, self.exit),
                       MessageHandler(Filters.text, self.exit)],
        )

        dispatcher.add_handler(conv_search, CONV_SEARCH)
