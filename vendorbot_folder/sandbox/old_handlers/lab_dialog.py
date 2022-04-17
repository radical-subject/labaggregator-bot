import logging

from bson import ObjectId
from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import CallbackContext, CommandHandler, CallbackQueryHandler, MessageHandler, Filters

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.decorators import log_errors
from modules.db import dbmodel


logger = logging.getLogger(__name__)


class LabDialog(Handlers):
    def __init__(self, db_instances):
        super().__init__(db_instances)
        
    def dialogue(self, update: Update, context: CallbackContext) -> None:
        current_state = context.user_data.get('state')
        input_text = update.message.text
        if current_state is None:
            return

        user_id = update.message.from_user.id
        chat_id = update.message.chat_id
        context.chat_data["user_id"] = user_id
        context.chat_data["chat_id"] = chat_id

        if current_state.startswith('LAB'):
            step = current_state.split(':')[1]
            if step.lower() == 'new':
                update.message.reply_text(f'New lab name is {input_text} \n\r Please enter organization name')
                context.user_data['state'] = 'LAB:NEW1'
                context.user_data['new_lab_name'] = input_text
                return self.LABS
            elif step.lower() == 'new1':
                update.message.reply_text(f'New lab organization is {input_text} \n\r Please enter contact')
                context.user_data['state'] = 'LAB:NEW2'
                context.user_data['new_lab_organization'] = input_text
                return self.LABS
            elif step.lower() == 'new2':
                update.message.reply_text(f'New lab contact  is {input_text}')
                context.user_data['state'] = ''
                context.user_data['new_lab_contact'] = input_text
                lab_dict = {
                    "name": context.user_data['new_lab_name'],
                    "organization": context.user_data['new_lab_organization'],
                    "contacts": context.user_data['new_lab_contact'],
                    "personal_list": [update.message.from_user.id]
                }
                context.user_data['current_lab'] = self.create_lab(lab_dict)
                context.user_data.pop('new_lab_name', None)
                context.user_data.pop('new_lab_organization', None)
                context.user_data.pop('new_lab_contact', None)
                return self.LABS
        elif current_state.startswith('SEARCH'):
            step = current_state.split(':')[1]
            if step.lower() == 'new':
                update.message.reply_text(
                    """
                    ⏳ 💅🏻\nSearching... please wait several seconds.\nИщу в базе... пожалуйста ожидайте. 
                    """)
                context.user_data['state'] = 'SEARCH:WAIT'

                try:
                    # RESOLVING DESCRIPTORS WITH CIRPY MODULE
                    SMILES = resolver.get_SMILES(input_text)

                    context.chat_data['SMILES'] = SMILES
                    context.chat_data['user_original_request_query'] = input_text
                    drawingPNGmodule.generate_png(SMILES, user_id)
                    img = open('./tmp/PNG/{}.png'.format(user_id), 'rb')
                    context.bot.send_photo(chat_id=chat_id, photo=img)
                    img.close()
                    button_list = [
                        [
                            InlineKeyboardButton("Yup", callback_data=str('YUP')),
                            InlineKeyboardButton("Nope", callback_data=str('NOPE'))
                        ]
                    ]
                    reply_markup = InlineKeyboardMarkup(button_list)
                    update.message.reply_text(
                        "😳 Did you mean this compound?\nТы искал(а) это соединение?",
                        reply_markup=reply_markup
                    )
                    return
                    # IN CASE IF NOTHING FOUND BY CIRPY OR QUERY IS STUPID
                except ValueError:  # in case when resolvingModule.resolve_query(request_query) = None, drawingSVGmodule.generate_picture_formula(SMILES) returns ValueError
                    update.message.reply_text(
                        """
                        Nothing found! Try to reformulate your request, baka.. Type 'exit' to quit search or send me another <CAS> or <chemical name> of desired reagent:
                    Ничего не нашлось! 🤷🏻‍♀️ Переформулируй запрос, или напиши мне "выход" чтобы закончить поиск:
                        """)
                return

    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(MessageHandler(Filters.all, self.dialogue))

    @log_errors
    def create_lab(self, lab_dict:dict):
        collection = 'laboratories'
        record = dbmodel.add_records(self.timerbot_db_client, self.db_instances["timerbot_db"], collection, lab_dict)
        lab = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], collection, {"_id": record.inserted_id})
        return lab
