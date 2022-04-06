import logging
import io
from telegram import Update, InlineKeyboardMarkup, ParseMode
from telegram.ext import (CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler)

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.decorators import log_errors
from modules.ourbot.handlers.helpers import get_txt_content

from modules.db import dbmodel, dbschema

logger = logging.getLogger(__name__)

UPLOAD_STATE = range(1)


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
        self.bot = bot
        self.collection = "users_collection"

    @log_errors
    def manage_entrypoint(self, update: Update, context: CallbackContext):
        update.message.reply_text('Отправьте мне .txt файл со списком CAS-номеров столбиком, '
                                  'следующего формата:\n\n<b>12411-12-3</b>\n<b>45646-23-2</b>\netc.\n\n'
                                  'Send cas list in .txt format.',
                                  parse_mode=ParseMode.HTML)

        return UPLOAD_STATE

    @log_errors
    def getting_file(self, update: Update, context: CallbackContext):

        # retrieving data from user message
        # ищем запись относящуюся к пользователю
        user_id = update.message.from_user.id
        mongo_query = {"user_id": user_id}
        user_info = update.message.from_user

        update.message.reply_text(f'Ожидайте: список обрабатывается.\nBe patient; it may take a while...')
        # Достаем из базы весь объект пользователя с реагентами
        # Если такого пользователя нет - функция на лету его создает и не плюется ошибками

        #TODO replace: user_reagents_object = dbschema.UserReagents(dbmodel.get_user(user_id))
        user_reagents_object = dbmodel.get_user_reagents_object(self.vendorbot_db_client, self.db_instances["vendorbot_db"], self.collection, mongo_query, user_info)

        current_state = context.user_data.get('state')

        CAS_list = get_txt_content(update, context)

        # оставляю возможность хардкодить вручную контакт, прописывая первую строку импортируемого файла руками:
        # в формате reagents_contact:+79265776746
        if CAS_list and CAS_list[0] and CAS_list[0].startswith("reagents_contact:"):
            contact = CAS_list[0].split(":")[1]
        else:
            # TODO удалить!!
            # его может не быть если он скрыт от ботов. это в настройках телеграма.
            # мы обрабатываем txt файл. НЕ нужно делать лишнюю работу.
            # я сделал возможным и без заполненного contact выводить:  contact ИЛИ username ИЛИ chat_id
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
        dbmodel.update_record(self.vendorbot_db_client, self.db_instances["vendorbot_db"], self.collection, mongo_query, data)

        # update.message.reply_text(f"{user_reagents_object.shared_reagents()[0]}")
        sent_message = update.message.reply_text(f'''file was successfully parsed and uploaded.
<b>import results</b>:
Строк в вашем списке: <b>{import_stats["input_lines_number"]}</b>
Правильных CAS-номеров: <b>{len(import_stats["valid_CAS_numbers"])}</b>
Опечатка в CAS: <b>{import_stats["failed_CAS_check_number"]}</b>
Не найдено SMILES для: <b>{import_stats["SMILES_not_found"]}</b> позиций
Найдено SMILES для: <b>{import_stats["SMILES_found"]}</b> реагентов
Прекурсоров найдено и вычеркнуто: <b>{import_stats["blacklist_filter_result"]}</b>

Итого: импортировано в базу <b>{import_stats["total_reagents_imported"]}</b> реагентов.
В вашей базе сейчас: <b>{import_stats["total_reagents_count_in_DB"]}</b> реагентов.
        ''', parse_mode='HTML')

        update.message.reply_text(sent_message)
        return ConversationHandler.END

    @log_errors
    def exit(self, update: Update, context: CallbackContext):
        """
        handler for terminating all dialog sequences
        """
        try:
            query = update.callback_query
            if query == None:
                reply_markup = InlineKeyboardMarkup([])
                update.message.reply_text("Таймер прерван другой командой.\nВсе переменные состояний очищены, и вы в этом виноваты сами.", reply_markup=reply_markup)
            else:
                update.message.reply_text(f"""Выход из диалога /manage.""")
        except:
            update.message.reply_text(f"""Выход из диалога /manage.""")
            pass
        # now clear all cached data
        # clear assosiated with user data and custom context variables
        context.chat_data.clear()
        context.user_data.clear()

        return ConversationHandler.END
    
    # @log_errors
    # def write_to_db(self, update: Update, context: CallbackContext, comment, archived_status="False"):
    #     """
    #     функция пишет информацию о потраченных минутах таймера в базу данных. 
    #     """
    #     # ищем запись относящуюся к пользователю
    #     user_id = update.message.from_user.id
    #     mongo_query = {"user_id": user_id}

    #     # дебажим содержимое переменных
    #     print("==============================")
    #     logger.info(f"mongo_query = {mongo_query}")
    #     logger.info(f"context_user_data = {context.user_data}")
    #     print("==============================")

    #     # user_data может не иметь "current_category"
    #     # усли бы было обращение к user_data по ключу то при отсутствии ключа получалось бы KeyError. чтобы этого избежать используем .get
    #     current_category = str(context.user_data.get('current_category')) # gives None if there is no such key == when no caegories were yet created
    #     total_elapsed_minutes = context.user_data["total_elapsed_minutes"]

    #     # достаем ее из бд
    #     previous_records=dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query)
        
    #     # результат поиска может оказаться пустым
    #     if previous_records == [] or previous_records == None:
    #         timer_object = dbschema.TimerData(
    #             **{
    #                 "user_id": user_id
    #             }
    #         )
    #         # к минимальному объекту добавляем текущую запись таймера
    #         try:
    #             context.user_data['hashtag']
    #             timer_object.add_timerdata_entry(total_elapsed_minutes, comment, current_category, archived_status, **{'hashtag' : context.user_data['hashtag']})
    #         except KeyError: 
    #             timer_object.add_timerdata_entry(total_elapsed_minutes, comment, current_category, archived_status)
    #         data = timer_object.export()
    #         # записываем в базу новосозданный объект
    #         dbmodel.add_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, data)
    #     else:
    #         # если раньше у пользователя были записи то импортируем данные пользователя в объект таймера
    #         timer_object = dbschema.TimerData(
    #             **previous_records[0]
    #         )
    #         # манипулируем с объектом таймера добавляя в него запись 
    #         # update db timer_data record by upserting and replasing existing data dict with updated 
    #         try:
    #             context.user_data['hashtag']
    #             timer_object.add_timerdata_entry(total_elapsed_minutes, comment, current_category, archived_status, **{'hashtag' : context.user_data['hashtag']})
    #         except KeyError: 
    #             timer_object.add_timerdata_entry(total_elapsed_minutes, comment, current_category, archived_status)
    #         data = timer_object.export()
    #         # записываем экспортированный в словарь объект в базу
    #         dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)


    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('end', self.exit))

        self.manage_dialog = ConversationHandler(
            entry_points=[CommandHandler('manage', self.manage_entrypoint)],
            states={
                    UPLOAD_STATE: [
                        MessageHandler(Filters.attachment, self.getting_file, run_async=True)
                    ]
                },
            fallbacks=[
                    MessageHandler(Filters.regex('^Done$'), self.exit)
                ]
        )

        dispatcher.add_handler(self.manage_dialog, 1)
