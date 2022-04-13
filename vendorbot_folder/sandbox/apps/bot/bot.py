import os
import time

from telegram import Update, ParseMode, ReplyKeyboardRemove
from telegram.ext import Updater, CommandHandler, CallbackContext, ConversationHandler,\
    MessageHandler, Filters

from logger import log


def start(update: Update, context: CallbackContext) -> None:
    user = update.effective_user

    update.message.reply_text(
        f"""Привет""",
        parse_mode=ParseMode.HTML)

    #todo db.add_user(...)


def searching_job(context: CallbackContext) -> None:

    time.sleep(2)  # test

    job = context.job
    context.bot.send_message(job.context, text='Поиск закончен')


def search(update: Update, context: CallbackContext) -> None:
    chat_id = update.message.chat_id

    remove_job_if_exists(f'search{chat_id}', context)
    context.job_queue.run_once(searching_job, 0, context=chat_id, name=f'search{chat_id}')

    text = 'Ищем...'
    update.message.reply_text(text)


def parsing_job(context: CallbackContext) -> None:

    time.sleep(2) # test

    job = context.job
    context.bot.send_message(job.context, text='Beep!')


def parse(update: Update, context: CallbackContext) -> None:
    chat_id = update.message.chat_id

    remove_job_if_exists(f'parse{chat_id}', context)
    context.job_queue.run_custom(parsing_job, {}, context=chat_id, name=f'parse{chat_id}')

    text = 'Обработка списка началась!'
    update.message.reply_text(text)


def remove_job_if_exists(name: str, context: CallbackContext) -> bool:
    current_jobs = context.job_queue.get_jobs_by_name(name)
    if not current_jobs:
        return False
    for job in current_jobs:
        job.schedule_removal()
    return True


def conv(update: Update, context: CallbackContext) -> int:
    update.message.reply_text(f'start conv. state=1')
    return 1


def test_handler(update: Update, context: CallbackContext) -> int:
    update.message.reply_text(f'Thank you! {update.message.text}')
    return 1

def conv2(update: Update, context: CallbackContext) -> int:
    update.message.reply_text(f'start conv2. state=1')
    return 1

def test2_handler(update: Update, context: CallbackContext) -> int:
    update.message.reply_text(f'Thank you2222 {update.message.text}')
    return 1


def cancel(update: Update, context: CallbackContext) -> int:
    """Cancels and ends the conversation."""
    user = update.message.from_user
    update.message.reply_text(
        'Bye!', reply_markup=ReplyKeyboardRemove()
    )

    return ConversationHandler.END


def main() -> None:

    updater = Updater(os.getenv('BOT_TOKEN'))

    dispatcher = updater.dispatcher

    dispatcher.add_handler(CommandHandler("start", start))
    dispatcher.add_handler(CommandHandler("help", start))
    dispatcher.add_handler(CommandHandler("manage", parse))
    dispatcher.add_handler(CommandHandler("search", search))

    conv_handler = ConversationHandler(
        entry_points=[CommandHandler('conv', conv)],
        states={
            1: [MessageHandler(Filters.text & ~Filters.command, test_handler),  # ~Filters.command чтобы отбрасывать команды
                CommandHandler('cancel', cancel)],
        },
        fallbacks=[CommandHandler('cancel', cancel)],  # если никакой хэндлер не отработает
    )

    conv_handler2 = ConversationHandler(
        entry_points=[CommandHandler('conv2', conv2)],
        states={
            1: [MessageHandler(Filters.text & ~Filters.command, test2_handler),
                # ~Filters.command чтобы отбрасывать команды
                CommandHandler('cancel', cancel)],
        },
        fallbacks=[CommandHandler('cancel', cancel)],  # если никакой хэндлер не отработает
    )

    dispatcher.add_handler(conv_handler)
    dispatcher.add_handler(conv_handler2)

    updater.start_polling()

    updater.idle()


if __name__ == '__main__':
    main()
