import logging
import traceback

from telegram import Update, ParseMode
from telegram.ext import (CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler)

from modules.db.users import users_collection
from modules.db.dbschema import get_contact
from modules.db.unique_molecules import unique_molecules_collection
from modules.db.parser import parse_reagent_list

from . import run_async
from .helpers import bot_commands_text, CONV_MANAGE, UPLOAD_STATE, get_file, \
    get_contact_from_dataframe, file_to_dataframe, df_to_reagents

logger = logging.getLogger(__name__)


class Manage:
    """
    В этом диалоговом хендлере пользователь добавляет свой список реагентов. 
    Присылается текстовый файл с CAS номерами.
    """

    def manage(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        logger.info(f"manage({chat_id})")

        user_id = update.message.from_user.id
        user = users_collection.get_user(user_id)
        if not user:
            update.message.reply_text("Это почему тебя нет в БД?! Тыкни /start")
            return ConversationHandler.END

        if not get_contact(user):
            update.message.reply_text("У нас нет твоих контактов. Тыкни /start")
            return ConversationHandler.END

        update.message.reply_text("Отправьте мне .txt файл со списком CAS-номеров столбиком, "
                                  "следующего формата:\n"
                                  "\n"
                                  "<b>12411-12-3</b>\n"
                                  "<b>45646-23-2</b>\n"
                                  "etc.\n"
                                  "\n"
                                  "Либо таблицу Excel формата:\n"
                                  "<b>[ CAS | location | name | contact]</b>\n",
                                  parse_mode=ParseMode.HTML)
        return UPLOAD_STATE

    def getting_file(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        user_id = update.message.from_user.id
        logger.info(f"getting_file({chat_id})")

        try:
            # Достаем из базы весь объект пользователя с реагентами
            # Пользователь должен быть
            user = users_collection.get_user(user_id)
            if not user:
                update.message.reply_text("Это почему тебя нет в БД?! Тыкни /start")
                return ConversationHandler.END

            contact = get_contact(user)
            if not contact:
                update.message.reply_text("У нас нет твоих контактов. Тыкни /start")
                return ConversationHandler.END

            update.message.reply_text(f"Ожидайте: список обрабатывается.\nBe patient; it may take a while...")

            # читаем txt, csv, xlsx
            file, file_name = get_file(update, context)

            cas_df = file_to_dataframe(file, file_name)
            if cas_df is None:
                update.message.reply_text('Valid file extensions are [ .txt | .csv | .xlsx ]')
                return ConversationHandler.END

            file_contact, cas_df = get_contact_from_dataframe(cas_df)

            reagents = df_to_reagents(cas_df)
            if file_contact:
                reagents, text_report = parse_reagent_list(reagents, file_contact)
            else:
                reagents, text_report = parse_reagent_list(reagents, contact)

            unique_molecules_collection.register_reagents(reagents)

            users_collection.clear_reagents(user_id)
            users_collection.add_reagents(user_id, reagents)

            sent_message = text_report + f"""Итого: База реагентов перезаписана. \
                                             Содержит <b>{users_collection.reagents_count(user_id)}</b> реагентов."""
        except Exception as err:
            logger.error(traceback.format_exc())
            sent_message = "Ошибка обработки, лаборанты уже разбужены!"

        update.message.reply_text(sent_message, parse_mode=ParseMode.HTML)

        update.message.reply_text(bot_commands_text(chat_id))

        return ConversationHandler.END

    def exit(self, update: Update, context: CallbackContext):
        """
        handler for terminating all dialog sequences
        """
        chat_id = update.message.chat_id
        logger.info(f"manage.exit({chat_id})")
        update.message.reply_text(bot_commands_text(chat_id))

        return ConversationHandler.END

    def register_handler(self, dispatcher):

        conv_manage = ConversationHandler(
            entry_points=[CommandHandler("manage", self.manage)],
            states={
                    UPLOAD_STATE: [
                        MessageHandler(Filters.attachment, self.getting_file, run_async=run_async())
                    ]
                },
            fallbacks=[MessageHandler(Filters.command, self.exit),
                       MessageHandler(Filters.text, self.exit)],
        )

        dispatcher.add_handler(conv_manage, CONV_MANAGE)
