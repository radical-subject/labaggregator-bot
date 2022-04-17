import logging

from telegram import Update, InlineKeyboardMarkup, InlineKeyboardButton
from telegram.ext import Filters, CallbackContext, CommandHandler, MessageHandler, ConversationHandler

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.decorators import log_errors
from modules.db import dbmodel, dbschema

logger = logging.getLogger(__name__)


class CategoriesDialog(Handlers):
    def __init__(self,db_instances):
        super().__init__(db_instances)
        self.collection = "timer_data_collection"

    @log_errors
    def add_category(self, update: Update, context: CallbackContext):
        update.message.reply_text('Create category ....\n\r Enter category name')
        context.user_data['state'] = 'CATEGORY:NEW'
        return self.CATEGORIES

    @log_errors
    def choose_category(self, update, context):
        user_id = update.message.from_user.id
        result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, {"user_id": user_id})
        try:
            categories_list = result[0]['categories']
        except:
            categories_list = None

        # categories_list = None if categories_list is None else categories_list
        menu = self.build_categories_menu(categories_list)
        msg_text = '{:^50}'.format('<b>🛠  CHOOSING CATEGORIES:</b>')
        msg_text += '\n'
        msg_text += '{:=^50}'.format('=')
        msg_text += '\n'
        msg_text += '<i>Select one from below, \n' \
                    'or create a new category:</i>'
        reply_markup = InlineKeyboardMarkup(menu)


        update.message.reply_text(text=msg_text, reply_markup=reply_markup, parse_mode='HTML')
        return self.CATEGORIES

    @log_errors
    def build_categories_menu(self, categories_list: list):
        markup = []
        if categories_list:
            for item in categories_list:
                markup.append(
                    [InlineKeyboardButton(f'  🛠  {item}', callback_data=f'CATEGORY_ID:{item}', )]
                    )

        markup.append([InlineKeyboardButton("  ➕  Create CATEGORY", callback_data=f'CATEGORY_ID:NEW'),
                       InlineKeyboardButton(" ❌  Cancel and exit", callback_data=f'CATEGORY_ID:CANCEL')])
        
        markup.append(
            [InlineKeyboardButton(f'СБРОСИТЬ ВЫБРАННУЮ КАТЕГОРИЮ', callback_data=f'CATEGORY_ID:NONE', )]
            )
        return markup

    @log_errors
    def dialog(self, update: Update, context: CallbackContext):
        """
        из декларирования функции убрано
        -> None:

        функция оперирует состояниями, которые хранятся в chat_data["state"]
        """
        
        current_state = context.user_data.get('state')
        input_text = update.message.text
        if current_state is None:
            return
        user_id = update.message.from_user.id
        chat_id = update.message.chat_id
        context.chat_data["user_id"] = user_id
        context.chat_data["chat_id"] = chat_id

        if current_state.startswith('CATEGORY'):
            step = current_state.split(':')[1]
            if step.lower() == 'new':

                update.message.reply_text(f'CATEGORY NAME assigned: {input_text}\nCATEGORY active: {input_text}')
                context.user_data['state'] = ''

                context.user_data['current_category'] = self.create_category(user_id, input_text)
                return self.CATEGORIES

    @log_errors
    def register_handler(self, dispatcher):
        """
        добавлем хендлеры в диспетчер
        """
        dispatcher.add_handler(CommandHandler('choose_category', self.choose_category))
        dispatcher.add_handler(CommandHandler('add_category', self.add_category))
        dispatcher.add_handler(MessageHandler(Filters.text, self.dialog))

    @log_errors
    def create_category(self, user_id, input_text: str):
        """
        write in DB category data
        """
        # находим интересующую нас запись в бд по id юзера
        mongo_query = {"user_id": user_id}
        previous_records=dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query)
        
        # если ничего не нашлось по данному id пользователя в базе:
        if previous_records==[]:
            # создаем инстанс таймера из существующей записи в бд
            timer_object = dbschema.TimerData(
                **{
                    "user_id": user_id
                }
            )
            timer_object.categories = [f"{input_text}"]
            data = timer_object.export()
            record = dbmodel.add_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, data)
            
            # выдаем последнюю записанную категорию
            search_result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, {"_id": data["_id"]})
            last_category = search_result[0]['categories'][-1]
            return last_category
        
        else:
            # создаем инстанс таймера из существующей записи в бд
            timer_object = dbschema.TimerData(
                **previous_records[0]
            )
            # обновляем данные в инстансе таймера добавляя запись о времени
            # update db timer_data record by upserting and replasing existing data dict with updated
            if "categories" in timer_object.export().keys():
                # если у объекта таймера есть уже хотя бы одна категория
                if f"{input_text}" not in timer_object.categories:
                    timer_object.categories.append(f"{input_text}")
                    data = timer_object.export()
                    record = dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)
                    
                    # выдаем последнюю записанную категорию
                    search_result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, {"_id": data["_id"]})
                    last_category = search_result[0]['categories'][-1]
                    return last_category
                else:
                    # ничего не записываем в базу поскольку введенная категория дублирует существующую
                    return None
            else:
                # если у объекта таймера вообще нет никаких категорий - а значит нет и аттрибута .categories
                timer_object.categories = [f"{input_text}"]
                data = timer_object.export()
                record = dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)

                # выдаем послднюю записанную категорию
                search_result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, {"_id": data["_id"]})
                last_category = search_result[0]['categories'][-1]
                return last_category
