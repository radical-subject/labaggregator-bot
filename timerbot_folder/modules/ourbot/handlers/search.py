import logging

from faker import Faker
from telegram import Update, InlineKeyboardMarkup, InlineKeyboardButton
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors
from modules.db import dbmodel

logger = logging.getLogger(__name__)


class Search(Handlers):
    def __init__(self, db_instances):
        super().__init__(db_instances)

    def add_to_wishlist(self, update: Update, context: CallbackContext):
        update.message.reply_text(text='add_to_wishlist')
        return

    def search(self, update: Update, context: CallbackContext):
        if context.user_data.get('current_lab') is None:
            update.message.reply_text(text='You dont selected laboratories '
                                           'At the first time use commnand /choose_lab ')
            return

        update.message.reply_text(
            "🙋🏻‍♀️ Enter query (name or CAS):\n\n🖋 Напиши название реагента (по-английски) или пришли интересующий CAS-номер:")
        context.user_data['state'] = 'SEARCH:NEW'
        return

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('add_to_wishlist', self.add_to_wishlist))
        dispatcher.add_handler(CommandHandler('search', self.search))
