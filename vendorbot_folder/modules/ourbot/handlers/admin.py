import logging, os
from io import BytesIO
from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler
from telegram.ext.dispatcher import run_async

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.decorators import log_errors, is_admin
from modules.ourbot.service import mongoDumpModule
from modules.db import dbconfig, dbmodel, rdkitdb, dbschema

logger = logging.getLogger(__name__)


class Admin(Handlers):
    def __init__(self, bot, db_instances):
        super().__init__(db_instances)
        # чтобы использовать модуль bot из пакета telegram здесь, 
        # нужно его передать при инициализации инстанса этого класса в OurBot.
        self.bot = bot
        self.collection = "users_collection"
    
    @log_errors
    @is_admin
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
    @is_admin
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
    @is_admin
    def digest(self, update: Update, context: CallbackContext):
        """
        produces digest
        """
        chat_id = update.message.chat_id

        update.message.reply_text(f'Ожидайте: список обрабатывается...')

        users = dbmodel.get_all_users() # dbmodel.iterate_over_collection_of_users(self.vendorbot_db_client, self.db_instances["vendorbot_db"], self.collection)
        # logger.info(f"len all_entries = {len(all_entries)}")

        digest = {}
        for user in users:
            for r in dbschema.get_shared_reagents(user):
                name = dbschema.reagent_name(r)
                contact = dbschema.reagent_contact(user, r)
                if name and name not in digest:
                    digest[name] = []

                if contact not in digest[name]:   # если у пользователя 2 реактива, то будет повторяться contact
                    digest[name].append(contact)

        logger.info(f"digest length = {len(digest.keys())}")
        update.message.reply_text(f'Всего {len(digest.keys())} CAS')

        if digest:
            cas_list = list(digest.keys())
            cas_list = sorted(cas_list)

            digest_cas_txt = '\n'.join(cas_list)

            f = BytesIO(bytes(digest_cas_txt, 'utf-8'))
            f.name = 'digest_cas.txt'
            f.seek(0)

            context.bot.send_document(chat_id, f)

            digest_txt = ''
            for cas in cas_list:
                digest_txt += f'{cas} : {", ".join(digest[cas])}\n'

            f = BytesIO(bytes(digest_txt, 'utf-8'))
            f.name = 'digest.txt'
            f.seek(0)

            context.bot.send_document(chat_id, f)

    @log_errors
    @is_admin
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

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('purge_handler', self.purge_handler))
        dispatcher.add_handler(CommandHandler('dump', self.dump, run_async=True))
        dispatcher.add_handler(CommandHandler('blacklist_update', self.update_rdkit_db_blacklist_handler, run_async=True))
        dispatcher.add_handler(CommandHandler('digest', self.digest, run_async=True))
