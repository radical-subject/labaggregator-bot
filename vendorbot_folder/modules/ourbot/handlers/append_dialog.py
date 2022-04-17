
from telegram import Update, ParseMode
from telegram.ext import (CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler)

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.helpers import get_txt_content

from modules.db.dbmodel import users_collection
from modules.db.dbschema import UserReagents, parse_cas_list
from modules.ourbot.handlers.helpers import bot_commands_text, CONV_APPEND, APPEND_STATE

from .manage_dialog import get_contact_from_cas_file

import logging
logger = logging.getLogger(__name__)


class Append(Handlers):
    """
    В этом диалоговом хендлере админ добавляет новые реагенты в свой список реагентов.
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

    def append(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        logger.info(f"append({chat_id})")
        update.message.reply_text("Отправьте мне .txt файл со списком CAS-номеров столбиком, "
                                  "следующего формата:\n\n<b>12411-12-3</b>\n<b>45646-23-2</b>\netc.\n\n"
                                  "Send cas list in .txt format.",
                                  parse_mode=ParseMode.HTML)

        return APPEND_STATE

    def append_getting_file(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        user_id = update.message.from_user.id
        user_info = update.message.from_user
        logger.info(f"getting_file({chat_id})")

        update.message.reply_text(f"Ожидайте: список обрабатывается.\nBe patient; it may take a while...")

        try:
            # Достаем из базы весь объект пользователя с реагентами
            # Пользователь должен быть
            user = users_collection.get_user(user_id)
            if not user:
                update.message.reply_text("Это почему тебя нет в БД?! Тыкни /start")
                return ConversationHandler.END

            user_reagents_object = UserReagents()

            cas_list = get_txt_content(update, context)

            contact, cas_list = get_contact_from_cas_file(cas_list, user_info)
            if not contact:
                update.message.reply_text("Добавьте первую строку 'reagents_contact:<телефон/почта>' "
                                          "или заполните свой username")
                return APPEND_STATE

            reagents, text_report = parse_cas_list(cas_list, contact)

            user_reagents_object.user_reagents += reagents  # добавляем

            data = user_reagents_object.export()

            users_collection.update_user(user_id, data)

            sent_message = text_report + f"""Итого: импортировано в базу <b>{len(reagents)}</b> реагентов.
                                             В вашей базе сейчас: <b>{len(user_reagents_object.user_reagents)}</b> реагентов."""

        except Exception as err:
            logger.error(err)
            sent_message = "Ошибка обработки, лаборанты уже разбужены!"

        update.message.reply_text(sent_message, parse_mode=ParseMode.HTML)

        update.message.reply_text(bot_commands_text(chat_id))

        return ConversationHandler.END

    def exit(self, update: Update, context: CallbackContext):
        """
        handler for terminating all dialog sequences
        """
        chat_id = update.message.chat_id
        logger.info(f"append.exit({chat_id})")
        update.message.reply_text(bot_commands_text(chat_id))

        context.chat_data.clear()
        context.user_data.clear()

        return ConversationHandler.END

    def register_handler(self, dispatcher):

        conv_append = ConversationHandler(
            entry_points=[CommandHandler("append", self.append)],
            states={
                APPEND_STATE: [
                    MessageHandler(Filters.attachment, self.append_getting_file, run_async=True)
                ]
            },
            fallbacks=[MessageHandler(Filters.command, self.exit),
                       MessageHandler(Filters.text, self.exit)],
        )

        dispatcher.add_handler(conv_append, CONV_APPEND)
