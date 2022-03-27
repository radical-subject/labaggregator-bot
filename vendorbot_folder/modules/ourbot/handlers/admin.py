import logging, os

from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler
from telegram.ext.dispatcher import run_async

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors, restricted
from modules.ourbot.service import mongoDumpModule
from modules.db import dbconfig, dbmodel, rdkitdb

logger = logging.getLogger(__name__)


class Admin(Handlers):
    def __init__(self, bot, db_instances):
        super().__init__(db_instances)
        # чтобы использовать модуль bot из пакета telegram здесь, 
        # нужно его передать при инициализации инстанса этого класса в OurBot.
        self.bot=bot
    
    @log_errors
    @restricted
    @run_async # deprecated way of async running # sometimes this may break! not now but configuration is dangerous overall
    def purge_handler(self, update: Update, context: CallbackContext):
        button_list = [
            [
                InlineKeyboardButton("✅ Yup", callback_data='ADMIN:YUP'),
                InlineKeyboardButton("🙅🏻‍♀️ Nope", callback_data='ADMIN:NOPE')
            ]
        ]
        reply_markup = InlineKeyboardMarkup(button_list)
        update.message.reply_text(
            "👩🏻‍🦰 Do you really want to purge the database?",
            reply_markup=reply_markup
        )
        return
    
    @log_errors
    @restricted
    def update_rdkit_db_blacklist_handler(self, update: Update, context: CallbackContext):
        """
        updates blacklist with srs/Narkotiki_test.sdf
        """
        reply = rdkitdb.update_rdkit_db_blacklist(self.blacklist_rdkit_db_client, self.db_instances["blacklist_rdkit_db"])
        update.message.reply_text(f"{reply} molecules successfully imported. nice!")
        reply = rdkitdb.update_blacklist_with_pandas (self.blacklist_rdkit_db_client, self.db_instances["blacklist_rdkit_db"])
        logging.info(f"{reply} molecules successfully imported with metadata in separate collection. nice!")
        return

        

    @log_errors
    def dump(self, update: Update, context: CallbackContext):
        """
        dumps whole db, archives it, and sends user .zip archive with data
        """
        chat_id = update.message.chat.id
        path = mongoDumpModule.dump_database(dbconfig.MONGO_INITDB_ROOT_USERNAME, dbconfig.MONGO_INITDB_ROOT_PASSWORD)[1]
        logging.info(f'{path}')
        files = os.listdir("./mongodumps")
        logging.info(f'{files}')
        # this bot cannot send more than 50 mb
        try:
            self.bot.sendDocument(chat_id=chat_id, document=open(f'{path}.zip', 'rb'), timeout=1000)
            result = 'Гена, помнишь ты просил меня принести тебе бекап базы данных, я пошел и принес, вот оно. Гена на.'
            update.message.reply_text(result)
        except:
            update.message.reply_text("что-то не так. скорее всего база данных слишком большая и тебе нужно наладить закачку на гуглодиск.")
        return

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('purge_handler', self.purge_handler))
        dispatcher.add_handler(CommandHandler('dump', self.dump))
        dispatcher.add_handler(CommandHandler('blacklist_update', self.update_rdkit_db_blacklist_handler))
