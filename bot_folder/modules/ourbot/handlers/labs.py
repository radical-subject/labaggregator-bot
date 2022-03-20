import logging

from faker import Faker
from telegram import Update, InlineKeyboardMarkup, InlineKeyboardButton
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors
from modules.db import dbmodel

logger = logging.getLogger(__name__)


class Labs(Handlers):
    def __init__(self,db_instances):
        super().__init__(db_instances)
        self.collection = "laboratories"

    def add_labs(self, update: Update, context: CallbackContext):
        update.message.reply_text('Lets go create labs ....\n\r Enter the name')
        context.user_data['state'] = 'LAB:NEW'
        return self.LABS

    def choose_lab(self, update, context):
        user_id = update.message.from_user.id
        result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, {"personal_list": user_id})
        result = 'not data' if result is None else result
        menu = self.build_labs_menu(result)
        msg_text = '{:^50}'.format('<b>🔬 Labs:</b>')
        msg_text += '\n'
        msg_text += '{:=^50}'.format('=')
        msg_text += '\n'
        msg_text += '<i>You can select one of the labs ' \
                    'or create new:</i>'
        reply_marlup = InlineKeyboardMarkup(menu)
        update.message.reply_text(text=msg_text, reply_markup=reply_marlup, parse_mode='HTML')
        return self.LABS

    def build_labs_menu(self, lst: list):
        markup = []

        for row in lst:
            markup.append(
                [InlineKeyboardButton(f'  🔬  {row.get("name")}', callback_data=f'LABS_ID:{str(row.get("_id"))}', )])

        markup.append([InlineKeyboardButton("  ➕  Create New Labs", callback_data=f'LABS_ID:NEW'),
                       InlineKeyboardButton(" ❌  Cancel Choose Labs", callback_data=f'LABS_ID:CANCEL')])
        return markup

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('choose_lab', self.choose_lab))
        dispatcher.add_handler(CommandHandler('add_lab', self.add_labs))
