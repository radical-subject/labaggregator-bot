import logging

from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import (Updater, CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler, CallbackQueryHandler)

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.service.decorators import log_errors
from modules.db import dbmodel, dbschema
from decimal import *

from modules.db.dbmodel import get_timerdata_object

logger = logging.getLogger(__name__)


class EditCategoriesDialog(Handlers):
    def __init__(self, bot, db_instances):
        super().__init__(db_instances)
        self.bot=bot
        self.collection = "timer_data_collection"
    
    @log_errors
    def edit_categories_entrypoint(self, update: Update, context: CallbackContext):
        '''
        launches dialog for editing categories:
        delete, archive, unarchive
        '''
        user_id = update.message.from_user.id
        result = dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, {"user_id": user_id})
        
        '''
        собираем лист со всеми категориями - сначала активные, потом - архивированные
        '''
        try:
            categories_list = result[0]['categories']
        except:
            categories_list = None

        try:
            categories_list += result[0]['archived_categories']
        except:
            pass

        '''
        теперь строим из них менюшку специальной функцией
        '''
        menu = self.build_categories_menu(categories_list)

        msg_text = '{:^50}'.format('<b>🛠   EDIT CATEGORIES:</b>')
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
    def build_categories_menu(self, categories_list: list):
        """
        вынесенный отдельно функционал построения менюшки с кнопками
        """
        markup = []
        if categories_list:
            for item in categories_list:
                markup.append(
                    [InlineKeyboardButton(f'  🛠  {item}', callback_data=f'CATEGORY_NAME:{item}', )]
                    )

        markup.append([InlineKeyboardButton(" ❌  Cancel and exit", callback_data=f'CATEGORY_SELECTION:CANCEL')])
        return markup


    @log_errors
    def select_action(self, update: Update, context: CallbackContext):
        """
        набор действий с категориями на выбор. 
        """
        current_state = context.user_data.get('state')
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        user_id = update.callback_query.from_user.id
        chat_id = update.callback_query.message.chat_id

        context.chat_data["user_id"] = user_id
        context.chat_data["chat_id"] = chat_id

        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message
        
        # достаем название категории из коллбека отдаваемого кнопкой
        category_name = query.data.split(":")[1]
        context.chat_data["category_name_manipulated"] = category_name

        # проверяем, в каком статусе находится категория - архивирована или нет? 
        result = dbmodel.get_records(
            self.timerbot_db_client,
            self.db_instances["timerbot_db"], 
            self.collection, 
            { 
                'user_id' : user_id, 
                'categories': { '$in': [f"{category_name}"] }
            }
        )

        """
        ниже альтернативный вариант поиска но менее надежный
        """
        # result = dbmodel.get_records(
        #     self.timerbot_db_client,
        #     self.db_instances["timerbot_db"], 
        #     self.collection, 
        #     { 
        #         "user_id" : user_id, 
        #         'timerdata': { '$elemMatch': {'category_name': { '$eq': category_name }, 'archived_status': { '$eq': "False" } } }
        #     }
        # )

        logger.info(result)

        """
        если категория не найдена в object.categories, значит она может быть только в archived_categories. 
        в зависимости от результата появляются кнопки архивировать либо разархивировать.
        """
        # logger.info(bool(result)) -- bool([]) == False !!
        # если не категория глобально не архивирована - result не будет None:
        if result:

            button_list = [
                [
                    InlineKeyboardButton("АРХИВИРОВАТЬ", callback_data=str('CATEGORY_ACTION:ARCHIVE')),
                    InlineKeyboardButton("УДАЛИТЬ", callback_data=str('CATEGORY_ACTION:DELETE_CONFIRM'))
                ]
            ]
            reply_markup = InlineKeyboardMarkup(button_list)

            self.bot.edit_message_text(
                text=f'выберите действие:',
                chat_id=sent_message.chat_id,
                message_id=sent_message.message_id,
                reply_markup=reply_markup
            )
        
        else:
            button_list = [
                [
                    InlineKeyboardButton("РАЗАРХИВИРОВАТЬ", callback_data=str('CATEGORY_ACTION:UNARCHIVE')),
                    InlineKeyboardButton("УДАЛИТЬ", callback_data=str('CATEGORY_ACTION:DELETE_CONFIRM'))
                ]
            ]
            reply_markup = InlineKeyboardMarkup(button_list)

            self.bot.edit_message_text(
                text=f'выберите действие:',
                chat_id=sent_message.chat_id,
                message_id=sent_message.message_id,
                reply_markup=reply_markup
            )


        return 2

    @log_errors
    def archive_category(self, update: Update, context: CallbackContext):
        """
        архивируем все записи связанные с этой категорией -
        - по сути меняем в них archived_status с False на True
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        category_name = context.chat_data["category_name_manipulated"]
        # ищем запись относящуюся к пользователю
        user_id = update.callback_query.from_user.id
        mongo_query = {"user_id": user_id}
        chat_id = update.callback_query.message.chat_id

        timerdata_object = get_timerdata_object(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, user_id)
        timerdata_object.mutate_archive_status(category_name)
        
        timerdata_object.categories.remove(f"{category_name}")
        try:
            timerdata_object.archived_categories.append(f"{category_name}")
        except:
            timerdata_object.archived_categories = [f"{category_name}"]

        data = timerdata_object.export()
        # записываем экспортированный в словарь объект в базу
        dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)


        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message
        # удаляем разметку кнопок
        reply_markup = InlineKeyboardMarkup([])

        self.bot.edit_message_text(
            text=f'категория архивирована.',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )

        return -1



    @log_errors
    def unarchive_category(self, update: Update, context: CallbackContext):
        """
        архивируем все записи связанные с этой категорией -
        - по сути меняем в них archived_status с False на True
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        category_name = context.chat_data["category_name_manipulated"]
        # ищем запись относящуюся к пользователю
        user_id = update.callback_query.from_user.id
        mongo_query = {"user_id": user_id}
        chat_id = update.callback_query.message.chat_id

        timerdata_object = get_timerdata_object(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, user_id)
        timerdata_object.mutate_archive_status(category_name)
        
        timerdata_object.archived_categories.remove(f"{category_name}")
        try:
            timerdata_object.categories.append(f"{category_name}")
        except:
            timerdata_object.categories = [f"{category_name}"]

        data = timerdata_object.export()
        # записываем экспортированный в словарь объект в базу
        dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)


        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message
        # удаляем разметку кнопок
        reply_markup = InlineKeyboardMarkup([])

        self.bot.edit_message_text(
            text=f'категория разархивирована.',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )

        return -1

    def delete_category_confirmation(self, update: Update, context: CallbackContext):
        button_list = [
            [
                InlineKeyboardButton("✅ Yup", callback_data='CATEGORY_ACTION:DELETE'),
                InlineKeyboardButton("🙅🏻‍♀️ Nope", callback_data='CATEGORY_ACTION:ABORT_DELETE')
            ]
        ]
        reply_markup = InlineKeyboardMarkup(button_list)
        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message

        self.bot.edit_message_text(
            text="👩🏻‍🦰 Do you really want to delete this category?\n(and all its related timer data)?",
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )

        return 3


    @log_errors
    def delete_category(self, update: Update, context: CallbackContext):
        """
        архивируем все записи связанные с этой категорией -
        - по сути меняем в них archived_status с False на True
        """
        query = update.callback_query
        # не забываем отвечать на query
        query.answer()

        category_name = context.chat_data["category_name_manipulated"]
        # ищем запись относящуюся к пользователю
        user_id = update.callback_query.from_user.id
        mongo_query = {"user_id": user_id}
        chat_id = update.callback_query.message.chat_id

        timerdata_object = get_timerdata_object(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, user_id)
        timerdata_object.delete_category(category_name)

        data = timerdata_object.export()
        # записываем экспортированный в словарь объект в базу
        dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)

        # достаем объект последнего отправленного ботом сообщения        
        sent_message = update.callback_query.message
        # удаляем разметку кнопок
        reply_markup = InlineKeyboardMarkup([])

        self.bot.edit_message_text(
            text=f'Категория удалена насовсем, включая все данные таймера о ней.',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )

        return -1



    @log_errors
    def exit(self, update: Update, context: CallbackContext):
        """
        handler for terminating all dialog sequences
        """
        try:
            query = update.callback_query
            if query == None:
                reply_markup = InlineKeyboardMarkup([])
                update.message.reply_text("Таймер абортирован другой командой.\nВсе переменные состояний очищены, и вы в этом виноваты сами.\nВыход из диалога таймера.\nДа поможет тебе святой Януарий!", reply_markup=reply_markup)
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

        self.edit_categories_dialog = ConversationHandler(
        entry_points=[CommandHandler('edit_categories', self.edit_categories_entrypoint)],
        states={
                1:[
                    CallbackQueryHandler(self.exit, pattern='^{}$'.format(str("CATEGORY_SELECTION:CANCEL"))),
                    CallbackQueryHandler(self.select_action, pattern='^{}'.format(str("CATEGORY_NAME")))
                ],
                2:[
                    CallbackQueryHandler(self.archive_category, pattern='^{}$'.format(str("CATEGORY_ACTION:ARCHIVE"))),
                    CallbackQueryHandler(self.unarchive_category, pattern='^{}$'.format(str("CATEGORY_ACTION:UNARCHIVE"))),
                    CallbackQueryHandler(self.delete_category_confirmation, pattern='^{}$'.format(str("CATEGORY_ACTION:DELETE_CONFIRM")))                                     
                ],
                3:[
                    CallbackQueryHandler(self.delete_category, pattern='^{}$'.format(str("CATEGORY_ACTION:DELETE"))),
                    CallbackQueryHandler(self.exit, pattern='^{}$'.format(str("CATEGORY_ACTION:ABORT_DELETE"))) 
                ]
            },
            fallbacks=[
                MessageHandler(Filters.regex('^Done$'), self.exit),
                MessageHandler(Filters.command, self.exit)
            ]
        )

        dispatcher.add_handler(self.edit_categories_dialog, 1)