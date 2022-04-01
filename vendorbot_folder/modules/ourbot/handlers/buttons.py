import logging, pymongo

from faker import Faker
from bson import ObjectId
from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler, CallbackQueryHandler

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors
from modules.db import dbmodel

logger = logging.getLogger(__name__)


class Buttons(Handlers):
    def __init__(self, db_instances):
        super().__init__(db_instances)

    @log_errors
    def button_handler(self, update: Update, context: CallbackContext) -> None:
        query = update.callback_query
        query.answer()
        faker = Faker()


        if query.data.startswith('LABS_ID'):
            """
            buttons for creation of laboratories are proceeded below
            """
            lab_id = query.data.split(':')[1]
            collection = 'laboratories'
            if lab_id == 'NEW':
                # update.message.reply_text('Lets go create labs ....\n\r Enter the name')
                context.user_data['state'] = 'LAB:NEW'
                query.edit_message_text(text=f"'Lets go create labs ....\n\r Enter the name'")
                return
            elif lab_id == 'CANCEL':
                context.user_data['state'] = 'LAB:CANCEL'
                query.edit_message_text(text=f"OK then.")
                return -1
            else:
                id = ObjectId(lab_id)
                result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], collection, {"_id": id})
                context.user_data['current_lab'] = result
                query.edit_message_text(text=f"Selected Lab: {result[0]}")
                return 

        elif query.data.startswith('CATEGORY_ID'): 
            """
            buttons for creation of categories are proceeded below
            """
            cat_name = query.data.split(':')[1]
            collection = 'timer_data_collection'
            if cat_name == 'NEW':
                context.user_data['state'] = 'CATEGORY:NEW'
                query.edit_message_text(text=f"Create category ...\n\r Enter category name")
                return
            elif cat_name == 'CANCEL':
                context.user_data['state'] = 'CATEGORY:CANCEL'
                query.edit_message_text(text=f"OK then.")
                return -1
            elif cat_name == 'NONE':
                context.user_data['current_category'] = None
                context.user_data['state'] = ''
                query.edit_message_text(text=f"Категория не выбрана = ведется запись дел второй важности")
                return
            else:
                context.user_data['current_category'] = cat_name
                query.edit_message_text(text=f"CATEGORY set active: {cat_name}")
                return


        elif query.data.startswith('ADMIN'):
            """
            buttons for database purging are proceeded below
            """
            command = query.data.split(':')[1]
            if command == 'YUP':
                try:
                    self.purge()
                    query.edit_message_text(text=f"База очищена.\nДа поможет тебе святой Франциск!")
                except pymongo.errors.OperationFailure:
                    query.edit_message_text(text=f"Client is not authorized on db to drop it.")

            if command == 'NOPE':
                query.edit_message_text(text=f"Блять нахуй я сюда пришёл... а, полотенце!")
        else:
            # query.edit_message_text(text=f"Something wrong")
            pass

    @log_errors
    def purge(self):
        """
        timerbot_client is not authorized on db to drop it.
        only root can. so, root_client and root instance is transferred as arguments to dbmodel.
        """
        dbmodel.purge(self.root_client, self.db_instances["vendorbot_db"])

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CallbackQueryHandler(self.button_handler))
        pass
