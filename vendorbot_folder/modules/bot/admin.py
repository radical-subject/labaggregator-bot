import os
import pymongo
import traceback
import logging

from io import BytesIO
from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler, CallbackQueryHandler

from modules.db.dbconfig import MONGO_VENDORBOT_DATABASE, root_client, MONGO_INITDB_ROOT_USERNAME, MONGO_INITDB_ROOT_PASSWORD

from modules.db.users import users_collection
from modules.db import dbschema
from modules.db.blacklist import blacklist_engine
from modules.db.mongodb import dump_database

from . import run_async
from .decorators import is_admin

logger = logging.getLogger(__name__)


class Admin:

    @is_admin
    def purge_handler(self, update: Update, context: CallbackContext):

        chat_id = update.message.chat_id
        logger.info(f"purge_handler({chat_id})")

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

    def button_handler(self, update: Update, context: CallbackContext):

        query = update.callback_query
        query.answer()

        logger.info(f"admin.button_handler(?) query data: {query.data}")

        if query.data == "ADMIN:YUP":
            try:
                self.purge()
                query.edit_message_text(text=f"База очищена.\nДа поможет тебе святой Франциск!")
            except pymongo.errors.OperationFailure:
                query.edit_message_text(text=f"Client is not authorized on db to drop it.")

        elif query.data == "ADMIN:NOPE":
            query.edit_message_text(text=f"Блять нахуй я сюда пришёл... а, полотенце!")

        else:
            pass  # сюда попадем из других callback, т.к. button_handler общий

    def purge(self):
        """
        vendorbot_db_client is not authorized on db to drop it.
        only root can. so, root_client and root instance is transferred as arguments to dbmodel.
        """
        root_client[MONGO_VENDORBOT_DATABASE].command("dropDatabase")
        logger.info(f"database {MONGO_VENDORBOT_DATABASE} dropped")

    @is_admin
    def blacklist_reload(self, update: Update, context: CallbackContext):
        """
        Загружает sdf файл со списком веществ в базу данных blacklist при помощи rdkit.
        updates blacklist with srs/Narkotiki_test.sdf
        """
        chat_id = update.message.chat_id
        logger.info(f"blacklist_reload({chat_id})")

        update.message.reply_text("Wait reloading...")
        reply = blacklist_engine.reload_rdkit()
        update.message.reply_text(f"{reply} molecules successfully imported. nice!")

        reply = blacklist_engine.reload_pandas()
        logger.info(f"{reply} molecules successfully imported with metadata in separate collection. nice!")

    @is_admin
    def digest(self, update: Update, context: CallbackContext):
        """
        Возвращает 2 файла:
        digest.txt – список shared компонентов в виде CAS и username пользователя
        digest_cas.txt - список всех shared компонентов в виде CAS
        """
        chat_id = update.message.chat_id
        logger.info(f"digest({chat_id})")

        update.message.reply_text(f'Ожидайте: список обрабатывается...')

        try:
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
            update.message.reply_text(f"Всего {len(digest.keys())} CAS")

            if digest:
                cas_list = list(digest.keys())
                cas_list = sorted(cas_list)

                digest_cas_txt = '\n'.join(cas_list)

                f = BytesIO(bytes(digest_cas_txt, "utf-8"))
                f.name = "digest_cas.txt"
                f.seek(0)

                context.bot.send_document(chat_id, f)

                digest_txt = '\n'.join(f'{cas} : {", ".join(digest[cas])}' for cas in cas_list)

                f = BytesIO(bytes(digest_txt, 'utf-8'))
                f.name = "digest.txt"
                f.seek(0)

                context.bot.send_document(chat_id, f)

        except Exception as err:
            logger.error(traceback.format_exc())
            update.message.reply_text(f'Ошибка формирования дайджеста...')

    @is_admin
    def dump(self, update: Update, context: CallbackContext):
        """
        dumps whole db, archives it, and sends user .zip archive with data
        """
        chat_id = update.message.chat.id
        logger.info(f'dump({chat_id})')

        response, path = dump_database(MONGO_INITDB_ROOT_USERNAME, MONGO_INITDB_ROOT_PASSWORD)
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
        dispatcher.add_handler(CommandHandler('dump', self.dump, run_async=run_async()))
        dispatcher.add_handler(CommandHandler('blacklist_reload', self.blacklist_reload, run_async=run_async()))
        dispatcher.add_handler(CommandHandler('digest', self.digest, run_async=run_async()))
        dispatcher.add_handler(CallbackQueryHandler(self.button_handler))
