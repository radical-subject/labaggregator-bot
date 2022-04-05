import logging, os

from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler
from telegram.ext.dispatcher import run_async

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors, restricted
from modules.ourbot.service import mongoDumpModule
from modules.db import dbconfig, dbmodel, rdkitdb, dbschema

logger = logging.getLogger(__name__)


class Admin(Handlers):
    def __init__(self, bot, db_instances):
        super().__init__(db_instances)
        # чтобы использовать модуль bot из пакета telegram здесь, 
        # нужно его передать при инициализации инстанса этого класса в OurBot.
        self.bot=bot
        self.collection = "users_collection"
    

    @log_errors
    @restricted
    # @run_async # deprecated way of async running # sometimes this may break! not now but configuration is dangerous overall
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
    @restricted
    def prepare_digest(self, update: Update, context: CallbackContext):
        """
        produces digest
        """
        
        # retrieving data from user message
        # ищем запись относящуюся к пользователю
        user_id = update.message.from_user.id
        mongo_query = {"user_id": user_id}
        user_info = update.message.from_user
        chat_id = update.message.chat.id

        update.message.reply_text(f'Ожидайте: список обрабатывается.\nBe patient; it may take a while...')
        # Достаем из базы весь объект пользователя с реагентами
        # Если такого пользователя нет - функция на лету его создает и не плюется ошибками
        all_entries = dbmodel.iterate_over_collection_of_users(self.vendorbot_db_client, self.db_instances["vendorbot_db"], self.collection)
        logger.info(f"len all_entries = {len(all_entries)}")
        didgest = None
        for entry in all_entries:
            user_reagents_object = dbschema.UserReagents(**entry)
            logger.info(len(user_reagents_object.user_reagents))
            # user_reagents_object = dbmodel.get_user_reagents_object(self.vendorbot_db_client, self.db_instances["vendorbot_db"], self.collection, mongo_query, user_info)
            # logger.info(user_reagents_object.export())
            didgest = user_reagents_object.get_digest_shared_reagents(didgest)
            logger.info(f"length = {len(didgest)}")
        
        logger.info(f"final length = {len(didgest)}")

        return 
        

    @log_errors
    @restricted
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
        dispatcher.add_handler(CommandHandler('dump', self.dump, run_async=True))
        dispatcher.add_handler(CommandHandler('blacklist_update', self.update_rdkit_db_blacklist_handler, run_async=True))
        dispatcher.add_handler(CommandHandler('didgest', self.prepare_digest, run_async=True))

