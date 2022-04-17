
from signal import SIGINT, SIGIOT, SIGPIPE

from telegram import Bot
from telegram.ext import Updater
from telegram.utils.request import Request

import logging
logger = logging.getLogger(__name__)
from modules.ourbot.handlers.initial import Inital
from modules.ourbot.handlers.admin import Admin
from modules.ourbot.handlers.search_dialog import Search
from modules.ourbot.handlers.manage_dialog import Manage
from modules.ourbot.handlers.append_dialog import Append


#def error_callback(update, context):
#    """
#    Унылый диалог - просто выводит сообщение и строку ошибки без трейслога
#    """
#    logger.warning('Update "%s" caused error "%s"', update, context.error)


class BotObject:

    def __init__(self, token: str, base_url: str = None, **db_instances):
        logger.info('Bot initialization... __init__ in BotObject...')
        self.token = token

        num_threads = 10
        request = Request(con_pool_size=num_threads + 4)
        self.bot = Bot(self.token, base_url, base_url, request=request)
        self.updater = Updater(bot=self.bot, workers=num_threads)

        self.dispatcher = self.updater.dispatcher
        self.db_instances = db_instances

        self.initial = Inital(self.bot, self.db_instances)
        self.admin = Admin(self.bot, self.db_instances)
        self.search_dialog = Search(self.bot, self.db_instances)
        self.manage_dialog = Manage(self.bot, self.db_instances)
        self.append_dialog = Append(self.bot, self.db_instances)

        logger.info('Bot initialization complete.')

    def start(self):
        logger.info('Starting bot...')
        self.update_dispatcher()
        self.updater.start_polling()
        self.updater.idle(stop_signals=(SIGINT, SIGIOT, SIGPIPE))  # TODO зачем другие сигналы?

    def update_dispatcher(self):
        #self.dispatcher.add_error_handler(error_callback)

        self.initial.register_handler(self.dispatcher)
        self.admin.register_handler(self.dispatcher)
        self.search_dialog.register_handler(self.dispatcher)
        self.manage_dialog.register_handler(self.dispatcher)
        self.append_dialog.register_handler(self.dispatcher)
