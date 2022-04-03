import logging
from telegram import (Bot, Update, InlineKeyboardMarkup, InlineKeyboardButton, ReplyKeyboardMarkup, KeyboardButton)
from telegram.ext import (Updater, CommandHandler, CallbackContext, ConversationHandler, InlineQueryHandler,
                          CallbackQueryHandler)
from telegram.utils.request import Request

from modules.ourbot.handlers.lab_dialog import LabDialog
from modules.ourbot.handlers.initial import Inital
from modules.ourbot.handlers.admin import Admin
from modules.ourbot.handlers.buttons import Buttons
from modules.ourbot.handlers.labs import Labs
from modules.ourbot.handlers.categories_dialog import CategoriesDialog
from modules.ourbot.handlers.edit_categories_dialog import EditCategoriesDialog
from modules.ourbot.handlers.search import Search
from modules.ourbot.handlers.wishlist import Wishlist
from modules.ourbot.handlers.manage_dialog import Manage

from modules.ourbot.handlers.timer_dialog import TimerDialog
from modules.ourbot.handlers.edit_entry_dialog import EditEntriesDialog

from modules.ourbot.service.decorators import log_errors
# from handlers.initial import register_initial_handler
logger = logging.getLogger(__name__)


class BotObject:

    def __init__(self, token: str, **db_instances):
        logger.info('Bot initialization... __init__ in BotObject...')
        self.token = token

        num_threads = 10
        request = Request(con_pool_size=num_threads + 4)
        self.bot = Bot(self.token, request=request)
        self.updater = Updater(bot=self.bot, workers=num_threads)

        # self.bot = Bot(self.token)
        # self.updater = Updater(self.token)
        self.dispatcher = self.updater.dispatcher
        self.db_instances=db_instances
        

        self.initial = Inital(self.bot, self.db_instances)
        self.admin = Admin(self.bot, self.db_instances)
        self.buttons = Buttons(self.db_instances)
        self.labs = Labs(self.db_instances)
        self.categories_dialog = CategoriesDialog(self.db_instances)
        self.edit_categories_dialog = EditCategoriesDialog(self.bot, self.db_instances)
        self.lab_dialog = LabDialog(self.db_instances)
        self.wishlist = Wishlist(self.db_instances)
        self.search = Search(self.db_instances)
        self.timer_dialog = TimerDialog(self.bot, self.db_instances)
        self.edit_entries_dialog = EditEntriesDialog(self.bot, self.db_instances)

        self.manage_dialog = Manage(self.bot, self.db_instances)
        


        logger.info('Bot initialization complete.')

    @log_errors
    def start(self):
        logger.info('Starting bot...')
        self.update_dispatcher()
        self.updater.start_polling()
        self.updater.idle(stop_signals=(2, 6, 13)) # stop_signals=(2, 6, 13)



    @log_errors
    def update_dispatcher(self):
        self.initial.register_handler(self.dispatcher)
        self.admin.register_handler(self.dispatcher)
        self.labs.register_handler(self.dispatcher)
        self.categories_dialog.register_handler(self.dispatcher)
        self.edit_categories_dialog.register_handler(self.dispatcher)
        self.wishlist.register_handler(self.dispatcher)
        self.search.register_handler(self.dispatcher)
        self.lab_dialog.register_handler(self.dispatcher)
        self.buttons.register_handler(self.dispatcher)
        self.timer_dialog.register_handler(self.dispatcher)
        self.edit_entries_dialog.register_handler(self.dispatcher)

        self.manage_dialog.register_handler(self.dispatcher)
 



        pass
