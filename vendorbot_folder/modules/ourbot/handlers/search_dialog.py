import os
import traceback
from telegram import Update, ReplyKeyboardMarkup, KeyboardButton, ReplyKeyboardRemove, ParseMode
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler, \
    MessageHandler, Filters

from modules.ourbot.handlers.helpers import CONV_SEARCH, SEARCH_STATE
from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.cas_to_smiles import pubchempy_smiles_resolve, cirpy_cas_resolve
from modules.ourbot.service.helpers import is_cas_number
import logging

from modules.db import dbschema
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
            def get_reagent_contact(users, text):
                ret = []
                if users:
                    for user in users:
                        reagents = dbschema.find_reagent(user, text)
                        for r in reagents:
                            contact = dbschema.reagent_contact(user, r)
                            if contact:
                                if contact not in contacts:
                                    ret.append(contact)
                return ret

            contacts = []

            if is_cas_number(text):
                update.message.reply_text("Ищем CAS в базе шеринга...")

                users = users_collection.get_users_by_cas(text)
                if users:
                    contacts.extend(get_reagent_contact(users, text))

            else:
                update.message.reply_text("Не похоже на CAS. Сейчас поищем по названию...")

                try:
                    users = users_collection.get_users_by_smiles(text)  # вдруг это smiles
                    if users:
                        contacts.extend(get_reagent_contact(users, text))

                    smiles = pubchempy_smiles_resolve(text)  # возвращает None на smiles
                    if smiles:
                        update.message.reply_text(f"Ищем по пользователям SMILES={smiles}")

                        users.extend(users_collection.get_users_by_smiles(smiles))
                        if users:
                            contacts.extend(get_reagent_contact(users, smiles))

                        cas = cirpy_cas_resolve(smiles)
                        if cas:
                            if isinstance(cas, str):  # TODO переписать cirpy_cas_resolve чтоыб всегда list возвращал
                                cas = [cas]

                            for c in cas:
                                update.message.reply_text(f"Ищем по пользователям CAS={c}")
                                users = users_collection.get_users_by_cas(c)
                                if users:
                                    contacts.extend(get_reagent_contact(users, c))
                        else:
                            update.message.reply_text(f"Не смогли определить CAS")

                except Exception as err:
                    logger.error(traceback.format_exc())
                    update.message.reply_text(f"Не смогли определить SMILES")
                    return SEARCH_STATE

            for user in users:
                reagents = dbschema.find_reagent(user, text)
                for r in reagents:
                    contact = dbschema.reagent_contact(user, r)
                    if contact:
                        if contact not in contacts:
                            contacts.append(contact)

            contacts = list(set(contacts))
            if contacts:
                update.message.reply_text(f'Реагентом могут поделиться эти контакты: {", ".join(contacts)}')
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
