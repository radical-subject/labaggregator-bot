import logging

from telegram import Update, InlineKeyboardButton, InlineKeyboardMarkup
from telegram.ext import (Updater, CommandHandler, MessageHandler, Filters, CallbackContext, ConversationHandler, CallbackQueryHandler)

from modules.ourbot.handlers.handlers import Handlers
from modules.ourbot.handlers.decorators import log_errors
from modules.db import dbmodel, dbschema
from modules.ourbot.service.timer import Timer, TimerError
from decimal import *

logger = logging.getLogger(__name__)


class TimerDialog(Handlers):
    def __init__(self, bot, db_instances):
        super().__init__(db_instances)
        self.bot=bot
        self.collection = "timer_data_collection"
    
    @log_errors
    def timer_entrypoint(self, update: Update, context: CallbackContext) -> None:
        """
        При начале диалога таймера таймер запускается сразу же.
        В зависимости от того, выбрана ли категория изначально, в конце работы таймера 
        формируется объект с информацией о:
        1) израсходованном времени (обязательно), 
        2) категории (опционально), 
        3) названии конкретного дела (обязательно) - заполняется после остановки таймера, 
        4) результатах (комментарий, опционально) - заполняется после остановки таймера,
        5) timestamp (обязательно). 
        """
        sent_message = update.message.reply_text(f'tic tac.')
        t=Timer()
        t.start()
        button_list = [
            [
                InlineKeyboardButton("||", callback_data=str('TIMER:PAUSE')),
                InlineKeyboardButton("|X|", callback_data=str('TIMER:STOP'))
            ]
        ]
        reply_markup = InlineKeyboardMarkup(button_list)
        self.bot.edit_message_text(
            text='Timer ON',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )
        context.user_data["timer_object"] = t

        return 1

    @log_errors
    def timer_pause(self, update: Update, context: CallbackContext) -> None:
        current_state = context.user_data.get('state')
        query = update.callback_query
        query.answer()
        if query.data.startswith('TIMER'):
            command = query.data.split(':')[1]
            assert command == "PAUSE", "wrong callback query data in timer_pause()"

        sent_message = update.callback_query.message
        t = context.user_data["timer_object"]
        seconds = t.stop()
        try:
            context.user_data["total_elapsed_seconds"] += seconds
        except:
            context.user_data["total_elapsed_seconds"] = seconds

        logger.info(context.user_data["total_elapsed_seconds"])
        context.user_data["total_elapsed_minutes"] = context.user_data["total_elapsed_seconds"]/60

        button_list = [
            [
                InlineKeyboardButton("|>", callback_data=str('TIMER:RESUME')),
                InlineKeyboardButton("|X|", callback_data=str('TIMER:STOP'))
            ]
        ]
        reply_markup = InlineKeyboardMarkup(button_list)
        self.bot.edit_message_text(
            text=f'Timer PAUSED.\ntime elapsed: {context.user_data["total_elapsed_minutes"]:.2f} min',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )
        context.user_data["timer_object"] = t

        return 2

    @log_errors
    def timer_continue(self, update: Update, context: CallbackContext) -> None:
        current_state = context.user_data.get('state')
        query = update.callback_query
        query.answer()
        if query.data.startswith('TIMER'):
            command = query.data.split(':')[1]
            assert command == "RESUME", "wrong callback query data in timer_pause()"

        sent_message = update.callback_query.message
        t = context.user_data["timer_object"]

        t.start()
        button_list = [
            [
                InlineKeyboardButton("||", callback_data=str('TIMER:PAUSE')),
                InlineKeyboardButton("|X|", callback_data=str('TIMER:STOP'))
            ]
        ]
        reply_markup = InlineKeyboardMarkup(button_list)
        self.bot.edit_message_text(
            text='Timer ON',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=reply_markup
        )
        context.user_data["timer_object"] = t

        return 1

    @log_errors
    def timer_stop(self, update: Update, context: CallbackContext) -> None:
        
        user_id = update.callback_query.from_user.id
        chat_id = update.callback_query.message.chat_id
        context.chat_data["user_id"] = user_id
        context.chat_data["chat_id"] = chat_id

        current_state = context.user_data.get('state')
        query = update.callback_query
        query.answer()
        if query.data.startswith('TIMER'):
            command = query.data.split(':')[1]
            assert command == "STOP", "wrong callback query data in timer_pause()"

        sent_message = update.callback_query.message
        t = context.user_data["timer_object"]
        try:
            seconds = t.stop()
            try:
                context.user_data["total_elapsed_seconds"] += seconds
            except:
                context.user_data["total_elapsed_seconds"] = seconds            
        except TimerError:
            logger.info("Timer was already stopped.")

        context.user_data["total_elapsed_minutes"] = context.user_data["total_elapsed_seconds"]/60
        
        logger.info(context.user_data["total_elapsed_seconds"])
        self.bot.edit_message_text(
            text=f'Timer STOPPED.\ntime elapsed: {context.user_data["total_elapsed_minutes"]:.2f} min.\n\n*Write down a comment:*',
            chat_id=sent_message.chat_id,
            message_id=sent_message.message_id,
            reply_markup=None,
            parse_mode='Markdown'
        )

        return 3
    
    @log_errors
    def send_data(self, update: Update, context: CallbackContext):
        input_text = update.message.text
        chat_id = update.message.chat_id
        comment = input_text
        user_id = update.message.from_user.id

        # self.bot.delete_message(chat_id=chat_id, message_id=update.message.message_id)
        # self.bot.delete_message(chat_id=chat_id, message_id=update.message.message_id-1)

        # write in DB context.user_data["total_elapsed_seconds"]
        current_category = str(context.user_data.get('current_category'))
        total_elapsed_minutes = context.user_data["total_elapsed_minutes"]

        # вызываем функцию для записи в бд
        self.write_to_db(update, context, comment)

        # ищем запись относящуюся к пользователю
        mongo_query = {"user_id": user_id}
        previous_records=dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query)
        timer_object = dbschema.TimerData(
            **previous_records[0]
        )

        try:
            time_netto_today = timer_object.get_netto_today()
            logger.info(f"time netto today = {timer_object.get_netto_today()} minutes")
        except:
            time_netto_today = 0
            
        # завершающее диалог сообщение пользователю об итогах таймера
        if current_category:
            logger.info(f"current_category = {current_category}")
            update.message.reply_text(f'Timer entry added successfully to database.\nCategory: {current_category}\nDone: *{comment}*\nTime spent: *{context.user_data["total_elapsed_minutes"]:.2f}* minutes. \n\ntime netto today = *{time_netto_today:.2f} minutes.*', parse_mode='Markdown')
        else:
            logger.info(f"current_category = {current_category}")
            update.message.reply_text(f'Timer entry added successfully to database.\nCategory was not chosen.\nDone: *{comment}*\nTime spent: *{context.user_data["total_elapsed_minutes"]:.2f}* minutes. \n\ntime netto today = *{time_netto_today:.2f} minutes.*', parse_mode='Markdown')

        # now clear all cached data
        # clear assosiated with user data and custom context variables
        context.chat_data.clear()
        context.user_data.clear()

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
            update.message.reply_text(f"""Выход из диалога.\nДа поможет тебе святой Франциск.""")
            pass
        # now clear all cached data
        # clear assosiated with user data and custom context variables
        context.chat_data.clear()
        context.user_data.clear()
        # equivalent of return ConversationHandler.END
        return -1
    
    @log_errors
    def write_to_db(self, update: Update, context: CallbackContext, comment, archived_status="False"):
        """
        функция пишет информацию о потраченных минутах таймера в базу данных. 
        """
        # ищем запись относящуюся к пользователю
        user_id = update.message.from_user.id
        mongo_query = {"user_id": user_id}

        # дебажим содержимое переменных
        print("==============================")
        logger.info(f"mongo_query = {mongo_query}")
        logger.info(f"context_user_data = {context.user_data}")
        print("==============================")

        # user_data может не иметь "current_category"
        # усли бы было обращение к user_data по ключу то при отсутствии ключа получалось бы KeyError. чтобы этого избежать используем .get
        current_category = str(context.user_data.get('current_category')) # gives None if there is no such key == when no caegories were yet created
        total_elapsed_minutes = context.user_data["total_elapsed_minutes"]

        # достаем ее из бд
        previous_records=dbmodel.get_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query)
        
        # результат поиска может оказаться пустым
        if previous_records == [] or previous_records == None:
            timer_object = dbschema.TimerData(
                **{
                    "user_id": user_id
                }
            )
            # к минимальному объекту добавляем текущую запись таймера
            try:
                context.user_data['hashtag']
                timer_object.add_timerdata_entry(total_elapsed_minutes, comment, current_category, archived_status, **{'hashtag' : context.user_data['hashtag']})
            except KeyError: 
                timer_object.add_timerdata_entry(total_elapsed_minutes, comment, current_category, archived_status)
            data = timer_object.export()
            # записываем в базу новосозданный объект
            dbmodel.add_records(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, data)
        else:
            # если раньше у пользователя были записи то импортируем данные пользователя в объект таймера
            timer_object = dbschema.TimerData(
                **previous_records[0]
            )
            # манипулируем с объектом таймера добавляя в него запись 
            # update db timer_data record by upserting and replasing existing data dict with updated 
            try:
                context.user_data['hashtag']
                timer_object.add_timerdata_entry(total_elapsed_minutes, comment, current_category, archived_status, **{'hashtag' : context.user_data['hashtag']})
            except KeyError: 
                timer_object.add_timerdata_entry(total_elapsed_minutes, comment, current_category, archived_status)
            data = timer_object.export()
            # записываем экспортированный в словарь объект в базу
            dbmodel.update_record(self.timerbot_db_client, self.db_instances["timerbot_db"], self.collection, mongo_query, data)

    @log_errors
    def set_hashtag(self, update: Update, context: CallbackContext):
        """
        """
        context.user_data['hashtag'] = update.message.text
        logger.info("DEBUG")
        logger.info(context.user_data['hashtag'])


    @log_errors
    def register_handler(self, dispatcher):
        dispatcher.add_handler(CommandHandler('end', self.exit))

        self.timer_dialog = ConversationHandler(
        entry_points=[CommandHandler('timer', self.timer_entrypoint)],
        states={
                1:[
                    CallbackQueryHandler(self.timer_pause, pattern='^{}$'.format(str("TIMER:PAUSE"))),
                    CallbackQueryHandler(self.timer_stop, pattern='^{}$'.format(str("TIMER:STOP"))),
                    MessageHandler(Filters.regex('^#[^ !@#$%^&*(),.?":{}|<>]*$'), self.set_hashtag)
                ],
                2:[
                    CallbackQueryHandler(self.timer_continue, pattern='^{}$'.format(str("TIMER:RESUME"))),
                    CallbackQueryHandler(self.timer_stop, pattern='^{}$'.format(str("TIMER:STOP"))),
                    MessageHandler(Filters.regex('^#[^ !@#$%^&*(),.?":{}|<>]*$'), self.set_hashtag)                
                ],
                3:[
                    # MessageHandler(Filters.command, self.exit),
                    MessageHandler(Filters.text, self.send_data)             
                ]
            },
            fallbacks=[
                MessageHandler(Filters.regex('^Done$'), self.exit)
                # MessageHandler(Filters.command, self.exit)
            ]
        )

        dispatcher.add_handler(self.timer_dialog, 1)