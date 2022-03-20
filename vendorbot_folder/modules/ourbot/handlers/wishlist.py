import logging

from telegram import Update, InlineKeyboardMarkup, InlineKeyboardButton
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors
from modules.db import dbmodel

logger = logging.getLogger(__name__)


class Wishlist(Handlers):
    def __init__(self, db_instances):
        super().__init__(db_instances)
        self.collection = "wishlists"

    def send(self, update: Update, context: CallbackContext):
        pass

    def list(self, update: Update, context: CallbackContext):
        pass

    def edit(self, update: Update, context: CallbackContext):
        pass

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('send', self.send))
        dispatcher.add_handler(CommandHandler('list', self.list))
        dispatcher.add_handler(CommandHandler('edit', self.edit))
