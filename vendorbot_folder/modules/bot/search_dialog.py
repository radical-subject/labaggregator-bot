import os
import traceback
import logging
from time import sleep

from telegram import Update, ReplyKeyboardMarkup, KeyboardButton, ReplyKeyboardRemove, error
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler, \
    MessageHandler, Filters

from . import run_async
from .helpers import CONV_SEARCH, SEARCH_STATE, is_admin_chat
from modules.chem.cas_to_smiles import what_reagent
from modules.chem import helpers
from modules.chem.pictures import create_similar_smiles_grid_picture, create_smiles_picture

from modules.db.dbschema import get_reagent_contacts, get_contact
from modules.db.users import users_collection
from modules.db.unique_molecules import *

logger = logging.getLogger(__name__)

CANCEL_SEARCH = 'Завершить поиск'
cancel_keyboard = [[KeyboardButton(CANCEL_SEARCH)]]

DBSIZE_OPEN_SEARCH = int(os.getenv('DBSIZE_OPEN_SEARCH', 10))


def send_photo(context, chat_id, file_path):
    context.bot.sendPhoto(chat_id=chat_id, photo=open(file_path, 'rb'), timeout=1000)


def unique_reagents(same_inchikey: str, same_smiles: str) -> List[Reagent]:
    reagents = users_collection.get_reagents_by_inchi(same_inchikey)
    # for smile_reagent in users_collection.get_reagents_by_smiles(same_smiles):
    #     # чтобы не повторялись (TODO не понятно reagent_id или inchi_key сравнивать)
    #     if not [r for r in reagents if r.inchikey_standard != smile_reagent.inchikey_standard]:
    #         reagents.append(smile_reagent)
    return reagents


def answer_user(user_id: int, reagents: List[Reagent]) -> str:
    """
    Отвечаем где находится реагенты:
    1. Если реагент у этого же пользователя, то вернем location (если есть)
    2. Если реагент у этого же пользователя, но заполнено поле contact, вернём его.
    3. Если реагент у другого пользователя, то вернём его контакт.
    :param user_id:
    :param reagents:
    :return:
    """
    logger.info([(r.contact, r.cas) for r in reagents])

    contacts = []
    locations = []
    for r in reagents:
        if r.user_id == user_id:
            if r.contact:
                contacts.append(r.contact)
            else:
                if r.location:
                    locations.append(r.location)
                else:
                    contacts.append('вас')
        else:
            contact = get_contact(users_collection.get_user(r.user_id))
            contacts.append(contact)

    locations = list(set(locations))
    contacts = list(set(contacts))

    reply_text = ''
    if locations:
        reply_text += f"Этот реагент есть у вас.\n" \
                      f"Попробуйте поискать его тут:\n{', '.join(locations)}."
        contacts.remove('вас')

    if contacts:
        if locations:
            reply_text += f"\nТакже этот реагент есть у {', '.join(contacts)}."
        else:
            reply_text += f"Этот реагент есть у {', '.join(contacts)}.\n"

    return reply_text


class Search:

    def search(self, update: Update, context: CallbackContext):
        """
        Старт ветки диалога "поиск"
        """
        chat_id = update.message.chat_id
        logger.info(f"search({chat_id})")

        user_id = update.message.from_user.id

        # Достаем из базы весь объект пользователя с реагентами
        # Пользователь должен быть
        user = users_collection.get_user(user_id)
        if not user:
            update.message.reply_text("Это почему тебя нет в БД?! Тыкни /start")
            return ConversationHandler.END

        count = len(users_collection.get_reagents(user_id))

        if count < DBSIZE_OPEN_SEARCH and not is_admin_chat(chat_id):
            update.message.reply_text(f"Чтобы разблокировать шеринг, вам необходимо загрузить "
                                      f"в базу не менее {DBSIZE_OPEN_SEARCH} ваших позиций. /manage")
            return ConversationHandler.END

        reply_markup = ReplyKeyboardMarkup(cancel_keyboard, resize_keyboard=True)
        update.message.reply_text("🙋🏻‍♀️ Enter query (name or CAS):\n\n"
                                  "🖋 Пришли интересующий **CAS-номер** или **название на английском языке**:",
                                  reply_markup=reply_markup)
        return SEARCH_STATE

    def search_cas(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        text = update.message.text
        user_id = update.message.from_user.id

        logger.info(f"search_cas({chat_id}): {text}")

        try:
            # ищем точное совпадение
            reagents = users_collection.get_reagents_by_text(text)
            if reagents:
                update.message.reply_text(answer_user(user_id, reagents))
            else:
                cas_list, smiles_list = what_reagent(text)
                if not cas_list and not smiles_list:
                    update.message.reply_text("Реагент не определён (ошибка в CAS?)")
                else:
                    ret = "Ищем по пользователям\n"
                    if cas_list:
                        ret += f"CAS: {', '.join(cas_list)}\n"
                    if smiles_list:
                        ret += f"SMILES: {', '.join(smiles_list)}"

                    update.message.reply_text(ret)

                    reagents_cas = users_collection.get_reagents_by_cas(cas_list)
                    reagents_smiles = users_collection.get_reagents_by_smiles(smiles_list)

                    reagents = list(set(reagents_cas + reagents_smiles))

                    if reagents:
                        update.message.reply_text(answer_user(user_id, reagents))
                    else:
                        for smiles in smiles_list:
                            filter_smile = helpers.filter_smiles_by_neutralize_atoms(smiles)
                            molecules = unique_molecules_collection.get_similar_molecules(filter_smile, 5)

                            if not molecules:
                                update.message.reply_text(f"Похожие на {smiles} не найдены")
                            else:  # что-то нашли
                                # similarity 0-1.0 схожесть где 1.0=100%
                                same_inchikey, same_smiles, similarity = molecules[0]
                                if similarity > 0.99:
                                    reagents_smiles = users_collection.get_reagents_by_smiles(same_smiles)
                                    reagents_inchikey = users_collection.get_reagents_by_inchi(same_inchikey)
                                    reagents = list(set(reagents_smiles + reagents_inchikey))

                                    if reagents:
                                        update.message.reply_text(f"Скорее всего вы ищете\n" \
                                                                  f"{reagents}" \
                                                                  f"(cхожесть с запросом {similarity*100}%):")
                                        update.message.reply_text(answer_user(user_id, reagents))

                                        path = create_smiles_picture(same_smiles)
                                        send_photo(context, chat_id, path)
                                    else:
                                        logger.error('Должны были найти, но не нашли: {same_smiles}')
                                        update.message.reply_text(f"Похожие на {same_smiles} не найдены")
                                else:
                                    update.message.reply_text(f"Точного совпадения нет.")

                                    for index, molecule in enumerate(molecules):
                                        same_inchikey, same_smiles, similarity = molecule

                                        reagents_smiles = users_collection.get_reagents_by_smiles(same_smiles)
                                        reagents_inchikey = users_collection.get_reagents_by_inchi(same_inchikey)
                                        reagents = list(set(reagents_smiles + reagents_inchikey))

                                        if reagents:
                                            update.message.reply_text(
                                                f"Совпадение с искомым реагентом "
                                                f"{similarity*100:.2f}% у:\n" \
                                                f"{', '.join([r.smiles for r in reagents])}")

                                            update.message.reply_text(answer_user(user_id, reagents))

                                        path = create_similar_smiles_grid_picture(smiles, molecules)
                                        send_photo(context, chat_id, path)

        except Exception as err:
            logger.error(traceback.format_exc())
            update.message.reply_text("Ошибка поиска. Похвастайтесь админу, что сломали бот.")

        logger.info(f"search_cas({chat_id}): {text} end")
        return SEARCH_STATE

    def exit(self, update: Update, context: CallbackContext) -> int:

        chat_id = update.message.chat_id
        logger.info(f'search.exit({chat_id})')

        update.message.reply_text("Поиск завершен",
                                  reply_markup=ReplyKeyboardRemove())

        return ConversationHandler.END

    def fallback(self, update: Update, context: CallbackContext) -> int:

        logger.error(f'{context}')

        return self.exit(update, context)

    def register_handler(self, dispatcher):

        conv_search = ConversationHandler(
            entry_points=[CommandHandler('search', self.search), ],
            states={
                SEARCH_STATE: [
                    MessageHandler(Filters.regex(CANCEL_SEARCH), self.exit),
                    MessageHandler(Filters.text & ~Filters.command, self.search_cas, run_async=run_async())
                ],
            },
            fallbacks=[MessageHandler(Filters.command, self.fallback),
                       MessageHandler(Filters.text, self.fallback)],
        )

        dispatcher.add_handler(conv_search, CONV_SEARCH)
