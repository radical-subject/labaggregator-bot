import logging
import traceback
from typing import List, Tuple, Optional

from telegram import Update, ParseMode
from telegram.ext import (CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler)

from modules.db.users import users_collection
from modules.db.dbschema import UserReagents, parse_cas_list, get_contact

from . import run_async
from .helpers import bot_commands_text, CONV_MANAGE, UPLOAD_STATE, get_file_content

import pandas as pd

logger = logging.getLogger(__name__)


def get_contact_from_first_row(df: pd.DataFrame) -> Tuple[Optional[str], pd.DataFrame]:
    """
    # оставляю возможность хардкодить вручную контакт для txt файлов,
    # прописывая первую строку импортируемого файла руками:
    # в формате reagents_contact:+79265776746
    :param df: содержимое файла
    :param user_info:
    :return:
    """
    if df['CAS'][0].startswith("reagents_contact:"):
        contact = df['CAS'][0]
        df = df.iloc[1:]
        return contact.strip("reagents_contact:"), df
    return None, df


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
                                  "следующего формата:\n\n<b>12411-12-3</b>\n<b>45646-23-2</b>\netc.\n\n"
                                  "Либо таблицу Excel формата:\n<b>[ CAS | location | name ]</b>\n\n"
                                  "Send cas list in .txt format.\n"
                                  "or Excel sheet with headers \n<b>[ CAS | location | name ]</b>",
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

            if not get_contact(user):
                update.message.reply_text("У нас нет твоих контактов. Тыкни /start")
                return ConversationHandler.END

            update.message.reply_text(f"Ожидайте: список обрабатывается.\nBe patient; it may take a while...")
            
            cas_tab = get_file_content(update, context) ### Remove after Pandas realisation

            if isinstance(cas_tab, str):
                update.message.reply_text(cas_tab)
                return ConversationHandler.END 
            #cas_tab = get_txt_content(update, context)
            

            contact, cas_list = get_contact_from_first_row(cas_tab)
            if contact:
                logger.info(f"getting_file({chat_id}): found contact in file {contact}")

            reagents, text_report = parse_cas_list(cas_list, contact)

            user_reagents_object = UserReagents(**user)

            user_reagents_object.user_reagents = reagents  # перезаписываем

            data = user_reagents_object.export()

            users_collection.update_user(user_id, data)

            sent_message = text_report + f"""Итого: База реагентов перезаписана. \
                                             Содержит <b>{len(user_reagents_object.user_reagents)}</b> реагентов."""

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
