
from typing import List, Tuple, Optional

from telegram import Update, ParseMode
from telegram.ext import (CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler)

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.helpers import get_txt_content

from modules.db.dbmodel import users_collection
from modules.db.dbschema import UserReagents, parse_cas_list, get_contact
from modules.ourbot.handlers.helpers import bot_commands_text, CONV_MANAGE, UPLOAD_STATE

import logging
import traceback
logger = logging.getLogger(__name__)


def get_contact_from_cas_file(cas_list: List[str]) -> Tuple[Optional[str], List[str]]:
    """
    # оставляю возможность хардкодить вручную контакт, прописывая первую строку импортируемого файла руками:
    # в формате reagents_contact:+79265776746
    :param cas_list: содержимое файла
    :param user_info:
    :return:
    """
    if cas_list and cas_list[0] and cas_list[0].startswith("reagents_contact:"):
        contact = cas_list[0].split(":")[1]
        cas_list = cas_list[1:]
        return contact, cas_list
    else:
        return None, cas_list


class Manage(Handlers):
    """
    В этом диалоговом хендлере пользователь добавляет свой список реагентов. 
    Присылается текстовый файл с CAS номерами.
    В качестве контакта обратной связи устанавливается контакт, который указан в шапке файла.
    Eсли не указан в шапке файла - то в качестве контакта используется username.

    Если у пользователя нет username - тогда все плохо. надо обработать этот момент. 
    """
    def __init__(self, bot, db_instances):
        """
        передаем коллекции, с которыми по умолчанию работает этот хендлер
        """
        super().__init__(db_instances)

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
                                  "Send cas list in .txt format.",
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

            cas_list = get_txt_content(update, context)

            contact, cas_list = get_contact_from_cas_file(cas_list)
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

        context.chat_data.clear()
        context.user_data.clear()
        return ConversationHandler.END

    def exit(self, update: Update, context: CallbackContext):
        """
        handler for terminating all dialog sequences
        """
        chat_id = update.message.chat_id
        logger.info(f"manage.exit({chat_id})")
        update.message.reply_text(bot_commands_text(chat_id))

        context.chat_data.clear()
        context.user_data.clear()
        return ConversationHandler.END

    def register_handler(self, dispatcher):

        conv_manage = ConversationHandler(
            entry_points=[CommandHandler("manage", self.manage)],
            states={
                    UPLOAD_STATE: [
                        MessageHandler(Filters.attachment, self.getting_file, run_async=True)
                    ]
                },
            fallbacks=[MessageHandler(Filters.command, self.exit),
                       MessageHandler(Filters.text, self.exit)],
        )

        dispatcher.add_handler(conv_manage, CONV_MANAGE)
