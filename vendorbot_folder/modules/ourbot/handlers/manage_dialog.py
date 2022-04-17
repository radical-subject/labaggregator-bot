
from telegram import Update, ParseMode
from telegram.ext import (CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler)

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.helpers import get_txt_content

from modules.db.dbmodel import users_collection
from modules.db import dbschema
from modules.ourbot.logger import logger
from modules.ourbot.handlers.helpers import bot_commands_text, CONV_MANAGE, UPLOAD_STATE


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
        logger.info(f'manage({chat_id})')
        update.message.reply_text('Отправьте мне .txt файл со списком CAS-номеров столбиком, '
                                  'следующего формата:\n\n<b>12411-12-3</b>\n<b>45646-23-2</b>\netc.\n\n'
                                  'Send cas list in .txt format.',
                                  parse_mode=ParseMode.HTML)

        return UPLOAD_STATE

    def getting_file(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        user_id = update.message.from_user.id
        user_info = update.message.from_user
        logger.info(f'getting_file({chat_id})')

        update.message.reply_text(f'Ожидайте: список обрабатывается.\nBe patient; it may take a while...')

        # Достаем из базы весь объект пользователя с реагентами
        # Пользователь должен быть
        user = users_collection.get_user(user_id)
        if not user:
            update.message.reply_text('Это почему тебя нет в БД?! Тыкни /start')
            return ConversationHandler.END

        user_reagents_object = dbschema.UserReagents()

        CAS_list = get_txt_content(update, context)

        # оставляю возможность хардкодить вручную контакт, прописывая первую строку импортируемого файла руками:
        # в формате reagents_contact:+79265776746
        if CAS_list and CAS_list[0] and CAS_list[0].startswith("reagents_contact:"):
            contact = CAS_list[0].split(":")[1]
        else:
            if not user_info.username:
                update.message.reply_text('Добавьте первую строку "reagents_contact:<телефон/почта>" '
                                          'или заполните свой username')
                return UPLOAD_STATE
            else:
                contact = f'@{user_info.username}'

        # импорт листа реагентов с фильтрациями
        import_stats = user_reagents_object.add_list_of_reagents(user_info.id, contact, self.blacklist_rdkit_db_client, self.db_instances["blacklist_rdkit_db"], CAS_list)
        # экспорт JSON - не работает с pymongo! нужен dict

        data = user_reagents_object.export()

        # записываем в базу объект
        users_collection.update_user(user_id, data)
        #mongo_query = {"user_id": user_id}
        #dbmodel.update_record(self.vendorbot_db_client, self.db_instances["vendorbot_db"], 'users_collection', mongo_query, data)

        sent_message = f'''file was successfully parsed and uploaded.
<b>import results</b>:
Строк в вашем списке: <b>{import_stats["input_lines_number"]}</b>
Правильных CAS-номеров: <b>{len(import_stats["valid_CAS_numbers"])}</b>
Опечатка в CAS: <b>{", ".join(import_stats["failed_CAS_check_number"])}</b>
Не найдено SMILES для: <b>{import_stats["SMILES_not_found"]}</b> позиций
Найдено SMILES для: <b>{import_stats["SMILES_found"]}</b> реагентов
Прекурсоров найдено и вычеркнуто: <b>{import_stats["blacklist_filter_result"]}</b>

Итого: импортировано в базу <b>{import_stats["total_reagents_imported"]}</b> реагентов.
В вашей базе сейчас: <b>{import_stats["total_reagents_count_in_DB"]}</b> реагентов.
        '''

        update.message.reply_text(sent_message, parse_mode=ParseMode.HTML)

        update.message.reply_text(bot_commands_text(chat_id))
        return ConversationHandler.END

    def exit(self, update: Update, context: CallbackContext):
        """
        handler for terminating all dialog sequences
        """
        chat_id = update.message.chat_id
        logger.info(f'manage.exit({chat_id})')
        update.message.reply_text(bot_commands_text(chat_id))

        context.chat_data.clear()
        context.user_data.clear()

        return ConversationHandler.END

    def register_handler(self, dispatcher):

        conv_manage = ConversationHandler(
            entry_points=[CommandHandler('manage', self.manage)],
            states={
                    UPLOAD_STATE: [
                        MessageHandler(Filters.attachment, self.getting_file, run_async=True)
                    ]
                },
            fallbacks=[MessageHandler(Filters.command, self.exit),
                       MessageHandler(Filters.text, self.exit)],
        )

        dispatcher.add_handler(conv_manage, CONV_MANAGE)
