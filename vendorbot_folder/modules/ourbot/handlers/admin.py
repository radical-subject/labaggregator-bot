import os
import traceback
from io import BytesIO
from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler

from modules.ourbot.logger import logger
from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.decorators import is_admin
from modules.ourbot.service import mongoDumpModule
from modules.db.dbmodel import users_collection
from modules.db import dbconfig, rdkitdb, dbschema


class Admin(Handlers):
    def __init__(self, bot, db_instances):
        super().__init__(db_instances)
        # чтобы использовать модуль bot из пакета telegram здесь, 
        # нужно его передать при инициализации инстанса этого класса в OurBot.
        self.bot = bot

    @is_admin
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

    @is_admin
    def update_rdkit_db_blacklist_handler(self, update: Update, context: CallbackContext):
        """
        Загружает sdf файл со списком веществ в базу данных blacklist при помощи rdkit.
        updates blacklist with srs/Narkotiki_test.sdf
        """
        reply = rdkitdb.update_rdkit_db_blacklist(self.blacklist_rdkit_db_client, self.db_instances["blacklist_rdkit_db"])
        update.message.reply_text(f"{reply} molecules successfully imported. nice!")
        reply = rdkitdb.update_blacklist_with_pandas (self.blacklist_rdkit_db_client, self.db_instances["blacklist_rdkit_db"])
        logger.info(f"{reply} molecules successfully imported with metadata in separate collection. nice!")

    @is_admin
    def digest(self, update: Update, context: CallbackContext):
        """
        Возвращает 2 файла:
        digest.txt – список shared компонентов в виде CAS и username пользователя
        digest_cas.txt - список всех shared компонентов в виде CAS
        """
        chat_id = update.message.chat_id

        update.message.reply_text(f'Ожидайте: список обрабатывается...')

        users = users_collection.get_all_users()
        logger.info(f"users {len(users)}")

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

    @is_admin
    def dump(self, update: Update, context: CallbackContext):
        """
        dumps whole db, archives it, and sends user .zip archive with data
        """
        chat_id = update.message.chat.id
        path = mongoDumpModule.dump_database(dbconfig.MONGO_INITDB_ROOT_USERNAME, dbconfig.MONGO_INITDB_ROOT_PASSWORD)[1]
        logger.info(f'{path}')
        files = os.listdir("./mongodumps")
        logger.info(f'{files}')
        # this bot cannot send more than 50 mb
        try:
            context.bot.sendDocument(chat_id=chat_id, document=open(f'{path}.zip', 'rb'), timeout=1000)
            result = 'Гена, помнишь ты просил меня принести тебе бекап базы данных, я пошел и принес, вот оно. Гена на.'
            update.message.reply_text(result)
        except Exception as err:
            tb = traceback.format_exc()
            update.message.reply_text("что-то не так. скорее всего база данных слишком большая и тебе нужно наладить закачку на гуглодиск.")
            update.message.reply_text(f"ошибка: {tb}")

    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('purge_handler', self.purge_handler))
        dispatcher.add_handler(CommandHandler('dump', self.dump, run_async=True))
        dispatcher.add_handler(CommandHandler('blacklist_update', self.update_rdkit_db_blacklist_handler, run_async=True))
        dispatcher.add_handler(CommandHandler('digest', self.digest))
