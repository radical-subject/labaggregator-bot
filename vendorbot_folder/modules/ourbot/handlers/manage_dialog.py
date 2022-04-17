
import uuid
import time

from telegram import Update, ParseMode
from telegram.ext import (CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler)

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.helpers import get_txt_content

from modules.db.dbmodel import users_collection
from modules.db import dbschema
from modules.db.blacklist import blacklist_engine
import logging
logger = logging.getLogger(__name__)
from modules.ourbot.handlers.helpers import bot_commands_text, CONV_MANAGE, UPLOAD_STATE

from modules.ourbot.service.helpers import is_cas_number
from modules.ourbot.service.cas_to_smiles import banch_cas_to_smiles


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
        update.message.reply_text("Отправьте мне .txt файл со списком CAS-номеров столбиком, "
                                  "следующего формата:\n\n<b>12411-12-3</b>\n<b>45646-23-2</b>\netc.\n\n"
                                  "Send cas list in .txt format.",
                                  parse_mode=ParseMode.HTML)

        return UPLOAD_STATE

    def getting_file(self, update: Update, context: CallbackContext):
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

            user_reagents_object = dbschema.UserReagents()

            cas_list = get_txt_content(update, context)

            # оставляю возможность хардкодить вручную контакт, прописывая первую строку импортируемого файла руками:
            # в формате reagents_contact:+79265776746
            if cas_list and cas_list[0] and cas_list[0].startswith("reagents_contact:"):
                contact = cas_list[0].split(":")[1]
            else:
                if not user_info.username:
                    update.message.reply_text("Добавьте первую строку 'reagents_contact:<телефон/почта>' "
                                              "или заполните свой username")
                    return UPLOAD_STATE
                else:
                    contact = f"@{user_info.username}"

            # импорт листа реагентов с фильтрациями

            valid_cas_list = [r for r in cas_list if is_cas_number(r)]
            failed_cas = [r for r in cas_list if not is_cas_number(r)]
            cas_smiles_list = banch_cas_to_smiles(valid_cas_list)

            cas_smiles_whitelist = [cas_smile for cas_smile in cas_smiles_list if not blacklist_engine.is_similar(cas_smile[1])]

            reagents = []

            now = time.strftime("%d.%m.%Y %H:%M", time.localtime())

            for cas, smiles in cas_smiles_whitelist:
                reagents.append({
                    "reagent_internal_id": uuid.uuid4().hex,
                    "CAS": cas,
                    "SMILES": smiles,
                    "contact": contact,
                    "sharing_status": "shared",
                    "timestamp": now
                })

            user_reagents_object.user_reagents = reagents   # перезаписываем

            data = user_reagents_object.export()

            # записываем в базу объект
            users_collection.update_user(user_id, data)

            sent_message = f"""file was successfully parsed and uploaded.
    <b>import results</b>:
    Строк в вашем списке: <b>{len(cas_list)}</b>
    Правильных CAS-номеров: <b>{len(valid_cas_list)}</b>
    Опечатка в CAS: <b>{", ".join(failed_cas)}</b>
    Не найдено SMILES для: <b>{len(valid_cas_list) - len(cas_smiles_list)}</b> позиций
    Найдено SMILES для: <b>{len(cas_smiles_list)}</b> реагентов
    Прекурсоров найдено и вычеркнуто: <b>{len(cas_smiles_list) - len(cas_smiles_whitelist)}</b>
    
    Итого: импортировано в базу <b>{len(cas_smiles_whitelist)}</b> реагентов.
    В вашей базе сейчас: <b>{len(user_reagents_object.user_reagents)}</b> реагентов.
            """
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
