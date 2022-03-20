import logging

from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import (Updater, CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler, CallbackQueryHandler)

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors
from modules.db import dbmodel, dbschema
from decimal import *

from modules.db.dbmodel import get_timerdata_object

logger = logging.getLogger(__name__)


class EditEntriesDialog(Handlers):
    def __init__(self, bot, db_instances):
        super().__init__(db_instances)
        self.bot=bot
        self.collection = "timer_data_collection"
    
    @log_errors
    def edit_entries_entrypoint(self, update: Update, context: CallbackContext):
        '''
        launches dialog for editing timer entries:
        timestamp editing, 
        '''
        user_id = update.message.from_user.id
        query = {"user_id": user_id}

        result = get_timerdata_object(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, query, user_id)
        '''
        достаем лист со всеми записями таймера:
        '''
        try:
            timer_data = result.timerdata
        except:
            timer_data = []
        # делаем список названий кнопок
        entry_list=[i['timestamp'] + " -- " + i['comment'] + " -- " + i['timerdata_entry_id'] for i in timer_data]
        '''
        теперь строим из них менюшку специальной функцией
        '''
        menu = self.build_menu(entry_list)

        msg_text = '{:^50}'.format('<b>✏️ EDIT ENTRIES:</b>')
        msg_text += '\n'
        msg_text += '{:=^50}'.format('=')
        msg_text += '\n'
        msg_text += '<i>Select entry from below:</i>'
        reply_markup = InlineKeyboardMarkup(menu)

        '''
        посылаем менюшку
        '''
        update.message.reply_text(text=msg_text, reply_markup=reply_markup, parse_mode='HTML')

        return 1

    @log_errors
    def build_menu(self, entry_list: list):
        """
        вынесенный отдельно функционал построения менюшки с кнопками
        """
        markup = []
        if entry_list:
            for item in entry_list:
                timestamp = item.split(" -- ")[0]
                comment = item.split(" -- ")[1]
                button_text = " -- ".join([timestamp, comment])
                timerdata_entry_id = item.split(" -- ")[-1]
                markup.append(
                    [InlineKeyboardButton(f'  🗒    {button_text}', callback_data=f'ENTRY_NAME:{timerdata_entry_id}', )]
                    )

        markup.append([InlineKeyboardButton(" ❌    CANCEL AND EXIT", callback_data=f'ENTRY_SELECTION:CANCEL')])
        return markup


    @log_errors
    def select_action(self, update: Update, context: CallbackContext):
        """
        набор действий с категориями на выбор. 
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        user_id = update.callback_query.from_user.id
        chat_id = update.callback_query.message.chat_id

        context.chat_data["user_id"] = user_id
        context.chat_data["chat_id"] = chat_id

        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message
        
        # достаем timestamp из коллбека отдаваемого кнопкой
        timerdata_entry_id = query.data.split(":")[1]
        # и сохраняем его в переменную чата - она нам еще понадобится далее
        context.chat_data["timerdata_entry_id"] = timerdata_entry_id

        # достаем нужную запись таймера сразу из базы маскированием посредством projection query
        result = dbmodel.get_records(
            self.timerbot_db_client,
            self.db_instances["timerbot_db"], 
            self.collection, 
            { 
                "user_id" : user_id, 
                'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } }
            },
            { 'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } } }
        )
        logger.info(result)

        str_line = ''
        for i in (list(item.items()) for item in result[0]['timerdata']):
            str_line += "\n".join([': '.join(map(str, tuple_item)) for tuple_item in i])
            str_line += "\n\n" 
        entry = str_line.rstrip("\n")

        if result:

            button_list = [
                [InlineKeyboardButton(" 📅   EDIT TIMESTAMP", callback_data=str('SELECT_ACTION:TIMESTAMP'))],
                [InlineKeyboardButton(" 🕔   EDIT TIME ELAPSED", callback_data=str('SELECT_ACTION:TIME'))],
                [InlineKeyboardButton(" 💬   EDIT COMMENT", callback_data=str('SELECT_ACTION:COMMENT'))],
                [InlineKeyboardButton(" 🗂   EDIT CATEGORY", callback_data=str('SELECT_ACTION:CATEGORY'))],
                [InlineKeyboardButton(" ❌   EXIT", callback_data=f'SELECT_ACTION:CANCEL')]
            ]
            reply_markup = InlineKeyboardMarkup(button_list)

            self.bot.edit_message_text(
                text=f'{entry}\n\nвыберите действие:',
                chat_id=sent_message.chat_id,
                message_id=sent_message.message_id,
                reply_markup=reply_markup
            )
        

        return 2

    @log_errors
    def edit_timestamp(self, update: Update, context: CallbackContext):
        """
        todo HORS анализ живого языка
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message
        context.chat_data['sent_message'] = sent_message
       
        reply_markup = InlineKeyboardMarkup([])

        '''
        посылаем менюшку и запрос комментария
        '''
        self.bot.edit_message_text(
            text=f'Напиши новую метку времени.\nЛибо в формате *YYYYMMDDHHMMSS*\nЛибо в произвольном текстовом формате, а я попробую угадать какая дата и время имелась в виду:',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )

        return 6

    @log_errors
    def edit_time(self, update: Update, context: CallbackContext):
        """
        добавить проверку данных - число ли? 
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message
        context.chat_data['sent_message'] = sent_message
       
        reply_markup = InlineKeyboardMarkup([])

        '''
        посылаем менюшку и запрос комментария
        '''
        self.bot.edit_message_text(
            text=f'Напиши количество времени в минутах (числом):',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )

        return 5

    @log_errors
    def edit_timestamp2(self, update: Update, context: CallbackContext):
        """
        фторой хендлер редактирования метки времени к записи таймера
        """
        user_info = update.message.from_user
        user_id = user_info.id
        # достаем объект последнего отправленного ботом сообщения        
        sent_message = context.chat_data['sent_message']
        timerdata_entry_id = context.chat_data["timerdata_entry_id"]

        mongo_query = {"user_id": user_id}
        result = get_timerdata_object(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, user_id)
        
        # достаем timestamp из cообщения юзера
        user_input = update.message.text
        # logger.info(user_input)
        result.edit_timestamp(timerdata_entry_id, user_input)
        data = result.export()
        # записываем экспортированный в словарь объект в базу
        dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)


        '''
        далее идет целый фрагмент функции select_action
        '''

        # достаем нужную запись таймера сразу из базы маскированием посредством projection query
        result = dbmodel.get_records(
            self.timerbot_db_client,
            self.db_instances["timerbot_db"], 
            self.collection, 
            { 
                "user_id" : user_id, 
                'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } }
            },
            { 'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } } }
        )
        # logger.info(result[0]['timerdata'])

        str_line = ''
        for i in (list(item.items()) for item in result[0]['timerdata']):
            str_line += "\n".join([': '.join(map(str, tuple_item)) for tuple_item in i])
            str_line += "\n\n" 
        entry = str_line.rstrip("\n")

        if result:

            button_list = [
                [InlineKeyboardButton(" 📅   EDIT TIMESTAMP", callback_data=str('SELECT_ACTION:TIMESTAMP'))],
                [InlineKeyboardButton(" 🕔   EDIT TIME ELAPSED", callback_data=str('SELECT_ACTION:TIME'))],
                [InlineKeyboardButton(" 💬   EDIT COMMENT", callback_data=str('SELECT_ACTION:COMMENT'))],
                [InlineKeyboardButton(" 🗂   EDIT CATEGORY", callback_data=str('SELECT_ACTION:CATEGORY'))],
                [InlineKeyboardButton(" ❌   EXIT", callback_data=f'SELECT_ACTION:CANCEL')]
            ]
            reply_markup = InlineKeyboardMarkup(button_list)

            self.bot.edit_message_text(
                text=f'{entry}\n\nвыберите действие:',
                chat_id=sent_message.chat_id,
                message_id=sent_message.message_id,
                reply_markup=reply_markup
            )

        return 2



    @log_errors
    def edit_time2(self, update: Update, context: CallbackContext):
        """
        фторой хендлер редактирования комментария к записи таймера
        """
        user_info = update.message.from_user
        chat_id = update.message.chat.id
        user_id = user_info.id
        # достаем объект последнего отправленного ботом сообщения        
        sent_message = context.chat_data['sent_message']
        timerdata_entry_id = context.chat_data["timerdata_entry_id"]

        mongo_query = {"user_id": user_id}
        result = get_timerdata_object(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, user_id)
        
        # достаем количество минут из cообщения юзера
        user_input = update.message.text
        # logger.info(user_input)
        result.edit_elapsed_time(timerdata_entry_id, user_input)
        data = result.export()
        # записываем экспортированный в словарь объект в базу
        dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)



        '''
        далее идет целый фрагмент функции select_action
        '''

        # достаем нужную запись таймера сразу из базы маскированием посредством projection query
        result = dbmodel.get_records(
            self.timerbot_db_client,
            self.db_instances["timerbot_db"], 
            self.collection, 
            { 
                "user_id" : user_id, 
                'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } }
            },
            { 'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } } }
        )
        # logger.info(result[0]['timerdata'])

        str_line = ''
        for i in (list(item.items()) for item in result[0]['timerdata']):
            str_line += "\n".join([': '.join(map(str, tuple_item)) for tuple_item in i])
            str_line += "\n\n" 
        entry = str_line.rstrip("\n")

        if result:

            button_list = [
                [InlineKeyboardButton(" 📅   EDIT TIMESTAMP", callback_data=str('SELECT_ACTION:TIMESTAMP'))],
                [InlineKeyboardButton(" 🕔   EDIT TIME ELAPSED", callback_data=str('SELECT_ACTION:TIME'))],
                [InlineKeyboardButton(" 💬   EDIT COMMENT", callback_data=str('SELECT_ACTION:COMMENT'))],
                [InlineKeyboardButton(" 🗂   EDIT CATEGORY", callback_data=str('SELECT_ACTION:CATEGORY'))],
                [InlineKeyboardButton(" ❌   EXIT", callback_data=f'SELECT_ACTION:CANCEL')]
            ]
            reply_markup = InlineKeyboardMarkup(button_list)

            self.bot.edit_message_text(
                text=f'{entry}\n\nвыберите действие:',
                chat_id=sent_message.chat_id,
                message_id=sent_message.message_id,
                reply_markup=reply_markup
            )

        return 2



    @log_errors
    def edit_comment(self, update: Update, context: CallbackContext):
        """
        вход в ветку редактирования комментария к записи таймера
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message
        context.chat_data['sent_message'] = sent_message
        user_id = update.callback_query.from_user.id
        # result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, {"user_id": user_id})
       
        reply_markup = InlineKeyboardMarkup([])

        '''
        посылаем менюшку и запрос комментария
        '''
        self.bot.edit_message_text(
            text=f'Напиши новый комментарий для этой записи:',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )

        return 4

    @log_errors
    def edit_comment2(self, update: Update, context: CallbackContext):
        """
        фторой хендлер редактирования комментария к записи таймера
        """
        user_info = update.message.from_user
        chat_id = update.message.chat.id
        user_id = user_info.id
        # достаем объект последнего отправленного ботом сообщения        
        sent_message = context.chat_data['sent_message']

        mongo_query = {"user_id": user_id}
        result = get_timerdata_object(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, user_id)
        
        # достаем название категории из cообщения юзера
        new_comment = update.message.text
        timerdata_entry_id = context.chat_data["timerdata_entry_id"]
        # Наконец заменяем категорию
        result.edit_comment(timerdata_entry_id, new_comment)

        data = result.export()
        # записываем экспортированный в словарь объект в базу
        dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)

        '''
        далее идет целый фрагмент функции select_action
        '''

        # достаем нужную запись таймера сразу из базы маскированием посредством projection query
        result = dbmodel.get_records(
            self.timerbot_db_client,
            self.db_instances["timerbot_db"], 
            self.collection, 
            { 
                "user_id" : user_id, 
                'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } }
            },
            { 'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } } }
        )
        # logger.info(result[0]['timerdata'])

        str_line = ''
        for i in (list(item.items()) for item in result[0]['timerdata']):
            str_line += "\n".join([': '.join(map(str, tuple_item)) for tuple_item in i])
            str_line += "\n\n" 
        entry = str_line.rstrip("\n")

        if result:

            button_list = [
                [InlineKeyboardButton(" 📅   EDIT TIMESTAMP", callback_data=str('SELECT_ACTION:TIMESTAMP'))],
                [InlineKeyboardButton(" 🕔   EDIT TIME ELAPSED", callback_data=str('SELECT_ACTION:TIME'))],
                [InlineKeyboardButton(" 💬   EDIT COMMENT", callback_data=str('SELECT_ACTION:COMMENT'))],
                [InlineKeyboardButton(" 🗂   EDIT CATEGORY", callback_data=str('SELECT_ACTION:CATEGORY'))],
                [InlineKeyboardButton(" ❌   EXIT", callback_data=f'SELECT_ACTION:CANCEL')]
            ]
            reply_markup = InlineKeyboardMarkup(button_list)

            self.bot.edit_message_text(
                text=f'{entry}\n\nвыберите действие:',
                chat_id=sent_message.chat_id,
                message_id=sent_message.message_id,
                reply_markup=reply_markup
            )

        return 2


    @log_errors
    def edit_category(self, update: Update, context: CallbackContext):
        """
        вход в ветку редактирования категории записи таймера
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message

        user_id = update.callback_query.from_user.id
        result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, {"user_id": user_id})
        
        '''
        собираем лист со всеми категориями - только активные
        '''
        try:
            categories_list = result[0]['categories']
        except:
            categories_list = None

        '''
        теперь строим из них менюшку специальной функцией
        '''
        menu = self.build_categories_menu(categories_list)

        msg_text = '{:^50}'.format('<b>🛠   ВЫБЕРИ НОВУЮ КАТЕГОРИЮ ДЛЯ ЭТОЙ ЗАПИСИ:</b>')
        msg_text += '\n'
        msg_text += '{:=^50}'.format('=')
        msg_text += '\n'
        msg_text += '<i>Select entry from below:</i>'
        reply_markup = InlineKeyboardMarkup(menu)

        '''
        посылаем менюшку
        '''
        self.bot.edit_message_text(
            text=f'Выбери новую категорию для этой записи:',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )

        return 3

    @log_errors
    def build_categories_menu(self, categories_list: list):
        """
        вынесенный отдельно функционал построения менюшки с кнопками
        """
        markup = []
        if categories_list:
            for item in categories_list:
                markup.append(
                    [InlineKeyboardButton(f'  🛠  {item}', callback_data=f'CANDIDATE_FOR_CATEGORY_NAME:{item}', )]
                    )

        markup.append([InlineKeyboardButton(" ❌  CANCEL", callback_data=f'CATEGORY_SELECTION_FOR_ENTRY:CANCEL')])
        return markup


    @log_errors
    def edit_category2(self, update: Update, context: CallbackContext):
        """
        меняем категорию на указанную в callback query
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message

        user_id = update.callback_query.from_user.id
        mongo_query = {"user_id": user_id}
        result = get_timerdata_object(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, user_id)
        
        # достаем название категории из коллбека отдаваемого кнопкой
        category_name = query.data.split(":")[1]
        timerdata_entry_id = context.chat_data["timerdata_entry_id"]
        # Наконец заменяем категорию
        result.edit_category(timerdata_entry_id, category_name)

        data = result.export()
        # записываем экспортированный в словарь объект в базу
        dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)

        '''
        далее идет целый фрагмент функции select_action
        '''

        # достаем нужную запись таймера сразу из базы маскированием посредством projection query
        result = dbmodel.get_records(
            self.timerbot_db_client,
            self.db_instances["timerbot_db"], 
            self.collection, 
            { 
                "user_id" : user_id, 
                'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } }
            },
            { 'timerdata': { '$elemMatch': {'timerdata_entry_id': { '$eq': timerdata_entry_id } } } }
        )
        # logger.info(result[0]['timerdata'])

        str_line = ''
        for i in (list(item.items()) for item in result[0]['timerdata']):
            str_line += "\n".join([': '.join(map(str, tuple_item)) for tuple_item in i])
            str_line += "\n\n" 
        entry = str_line.rstrip("\n")

        if result:

            button_list = [
                [InlineKeyboardButton(" 📅   EDIT TIMESTAMP", callback_data=str('SELECT_ACTION:TIMESTAMP'))],
                [InlineKeyboardButton(" 🕔   EDIT TIME ELAPSED", callback_data=str('SELECT_ACTION:TIME'))],
                [InlineKeyboardButton(" 💬   EDIT COMMENT", callback_data=str('SELECT_ACTION:COMMENT'))],
                [InlineKeyboardButton(" 🗂   EDIT CATEGORY", callback_data=str('SELECT_ACTION:CATEGORY'))],
                [InlineKeyboardButton(" ❌   EXIT", callback_data=f'SELECT_ACTION:CANCEL')]
            ]
            reply_markup = InlineKeyboardMarkup(button_list)

            self.bot.edit_message_text(
                text=f'{entry}\n\nвыберите действие:',
                chat_id=sent_message.chat_id,
                message_id=sent_message.message_id,
                reply_markup=reply_markup
            )

        return 2


    @log_errors
    def exit(self, update: Update, context: CallbackContext):
        """
        handler for terminating all dialog sequences
        """
        try:
            query = update.callback_query
            if query == None:
                reply_markup = InlineKeyboardMarkup([])
                update.message.reply_text("Процесс абортирован другой командой.\nВсе переменные состояний очищены, и вы в этом виноваты сами.\nВыход из диалога редактирования записей.\nДа поможет тебе святой Януарий!", reply_markup=reply_markup)
            else:
                update.message.reply_text(f"""Выход из диалога.\nДа поможет тебе святой Антоний.""")
        except:
            sent_message = update.callback_query.message
            self.bot.edit_message_text(
                text=f'Выход из диалога.\nДа поможет тебе святой Серафим Арзамас-16ый.',
                chat_id=sent_message.chat_id,
                message_id=sent_message.message_id,
                reply_markup=None,
                parse_mode='Markdown'
            )
            pass
        # now clear all cached data
        # clear assosiated with user data and custom context variables
        context.chat_data.clear()
        context.user_data.clear()
        # equivalent of return ConversationHandler.END
        return -1


    @log_errors
    def register_handler(self, dispatcher):
        """
        регистрируем все описанные выше хендлеры
        """
        dispatcher.add_handler(CommandHandler('end', self.exit))

        self.edit_entries_dialog = ConversationHandler(
        entry_points=[CommandHandler('edit_entries', self.edit_entries_entrypoint)],
        states={
                1:[
                    CallbackQueryHandler(self.exit, pattern='^{}$'.format(str("ENTRY_SELECTION:CANCEL"))),
                    CallbackQueryHandler(self.select_action, pattern='^{}'.format(str("ENTRY_NAME")))
                ],
                2:[
                    CallbackQueryHandler(self.edit_timestamp, pattern='^{}$'.format(str("SELECT_ACTION:TIMESTAMP"))),
                    CallbackQueryHandler(self.edit_time, pattern='^{}$'.format(str("SELECT_ACTION:TIME"))),
                    CallbackQueryHandler(self.edit_comment, pattern='^{}$'.format(str("SELECT_ACTION:COMMENT"))),
                    CallbackQueryHandler(self.edit_category, pattern='^{}$'.format(str("SELECT_ACTION:CATEGORY"))),
                    CallbackQueryHandler(self.exit, pattern='^{}$'.format(str("SELECT_ACTION:CANCEL")))
                ],
                3:[
                    CallbackQueryHandler(self.edit_category2, pattern='^{}'.format(str("CANDIDATE_FOR_CATEGORY_NAME"))),
                    CallbackQueryHandler(self.exit, pattern='^{}$'.format(str("CATEGORY_SELECTION_FOR_ENTRY:CANCEL"))) 
                ],
                4:[
                    MessageHandler(Filters.text, self.edit_comment2)
                ],
                5:[
                    MessageHandler(Filters.regex('^\d*$'), self.edit_time2)
                ],
                6:[
                    MessageHandler(Filters.regex('^(\d{9}|\d{14})$'), self.edit_timestamp2)
                ]
                
            },
            fallbacks=[
                MessageHandler(Filters.regex('^Done$'), self.exit),
                MessageHandler(Filters.command, self.exit)
            ]
        )

        dispatcher.add_handler(self.edit_entries_dialog, 1)