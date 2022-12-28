import os
import traceback
import logging

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps

from telegram import Update, ReplyKeyboardMarkup, KeyboardButton, ReplyKeyboardRemove
from telegram.ext import CallbackContext, CommandHandler, ConversationHandler, \
    MessageHandler, Filters

from . import run_async
from .helpers import CONV_SEARCH, SEARCH_STATE
from modules.chem.cas_to_smiles import what_reagent

from modules.db.dbschema import get_reagent_contacts, get_contact
from modules.db.users import users_collection
from modules.db.unique_molecules import *

logger = logging.getLogger(__name__)

CANCEL_SEARCH = 'Завершить поиск'
cancel_keyboard = [[KeyboardButton(CANCEL_SEARCH)]]


DBSIZE_OPEN_SEARCH = int(os.getenv('DBSIZE_OPEN_SEARCH', 10))


class Search:

    def search(self, update: Update, context: CallbackContext):
        """
        Старт ветки диалога "поиск"
        """
        chat_id = update.message.chat_id
        logger.info(f"search({chat_id})")

        user_id = update.message.from_user.id
        count = len(users_collection.get_reagents(user_id))

        if count < DBSIZE_OPEN_SEARCH:  #  and not is_admin_chat(chat_id)
            update.message.reply_text(f"Чтобы разблокировать шеринг, вам необходимо загрузить "
                                      f"в базу не менее {DBSIZE_OPEN_SEARCH} ваших позиций. /manage")
            return ConversationHandler.END

        reply_markup = ReplyKeyboardMarkup(cancel_keyboard, resize_keyboard=True)
        update.message.reply_text("🙋🏻‍♀️ Enter query (name or CAS):\n\n"
                                  "🖋 Пришли интересующий CAS-номер:",
                                  reply_markup=reply_markup)

        return SEARCH_STATE

    def search_cas(self, update: Update, context: CallbackContext):
        chat_id = update.message.chat_id
        text = update.message.text

        logger.info(f"search_cas({chat_id}): {text}")

        try:
            contacts = []

            cas_list, smiles_list = what_reagent(text)

            if cas_list or smiles_list:
                text = "Ищем по пользователям:"
                if cas_list:
                    text += f"\nCAS: {', '.join(cas_list)}\n"
                if smiles_list:
                    text += f"\nSMILES: {', '.join(smiles_list)}"
                update.message.reply_text(text)

                for cas in cas_list:
                    contacts.extend(get_reagent_contacts(users_collection.get_users_by_cas(cas), cas))

                for smiles in smiles_list:
                    contacts.extend(get_reagent_contacts(users_collection.get_users_by_smiles(smiles), smiles))

                contacts = list(set(contacts))
                
                if contacts:
                    update.message.reply_text(f"Реагентом могут поделиться эти контакты: {', '.join(contacts)}")
                else:
                    """
                    ВНИМАНИЕ !!!!
                    ТУТ СТРАННОЕ МЕСТО, дебажить в первую очередь если поиск чудит
                    """
                    best_match_smiles = unique_molecules_collection.get_most_similar_reagent(smiles)
                    
                    if best_match_smiles != None:
                        # for reagent_id in best_match_smiles[0]:
                        inchi_key = best_match_smiles[0]
                        contacts += [get_contact(i) for i in users_collection.get_user_by_reagent_inchi_key(inchi_key)]
                        update.message.reply_text(f"Кажется, реагентом пока никто не готов поделиться, но найден похожий у {', '.join(contacts)}.\nсхожесть с запросом: {(best_match_smiles[2]*100):.2f}%\n{best_match_smiles[1:]}\nSimilarity Map Result:")
        
    
                        mol = Chem.MolFromSmiles(best_match_smiles[1])
                        refmol = Chem.MolFromSmiles(smiles)

                        fp = SimilarityMaps.GetAPFingerprint(mol, fpType='normal')
                        fp = SimilarityMaps.GetTTFingerprint(mol, fpType='normal')
                        fp = SimilarityMaps.GetMorganFingerprint(mol, fpType='bv')
                        fig, maxweight = SimilarityMaps.GetSimilarityMapForFingerprint(refmol, mol, SimilarityMaps.GetMorganFingerprint)
                        fig.savefig(f"/vendorbot_container/srs/pic/{inchi_key}.png", bbox_inches = "tight")
                        
                        path = "/vendorbot_container/srs/pic"
                        context.bot.sendPhoto(chat_id=chat_id, photo=open(f'{path}/{inchi_key}.png', 'rb'), timeout=1000)
                        # result = f'Similarity Map Result. \nсхожесть с запросом: {(best_match_smiles[1]*100):.2f}%'
                        # update.message.reply_text(result)

                    else: 
                        update.message.reply_text(f"Реагентом пока никто не готов поделиться")
            else:
                update.message.reply_text("Реагент не определен (ошибка в CAS?)")

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
            entry_points=[CommandHandler('search', self.search),],
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
