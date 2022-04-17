import os

from telegram import Update, ParseMode, ReplyKeyboardRemove
from telegram.ext import Updater, CommandHandler, CallbackContext, ConversationHandler,\
    MessageHandler, Filters


def start(update: Update, context: CallbackContext) -> None:
    user = update.effective_user

    update.message.reply_text(
        f"""Привет""",
        parse_mode=ParseMode.HTML)


def conv(update: Update, context: CallbackContext) -> int:
    update.message.reply_text(f'start conv. state=1')
    return 1


def test_handler(update: Update, context: CallbackContext) -> int:
    update.message.reply_text(f'test_handler {update.message.text}')
    return 1


def conv2(update: Update, context: CallbackContext) -> int:
    update.message.reply_text(f'start conv2. state=1')
    return 1


def test2_handler(update: Update, context: CallbackContext) -> int:
    update.message.reply_text(f'test2_handler {update.message.text}')
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

    conv_handler = ConversationHandler(
        entry_points=[CommandHandler('conv', conv)],
        states={
            1: [MessageHandler(Filters.text & ~Filters.command, test_handler),  # ~Filters.command чтобы отбрасывать команды
                CommandHandler('cancel', cancel)],
        },
        fallbacks=[MessageHandler(Filters.command, cancel)],  # если никакой хэндлер не отработает
    )

    conv_handler2 = ConversationHandler(
        entry_points=[CommandHandler('conv2', conv2)],
        states={
            1: [MessageHandler(Filters.text & ~Filters.command, test2_handler),
                CommandHandler('cancel', cancel)],
        },
        fallbacks=[MessageHandler(Filters.command, cancel)],  # если никакой хэндлер не отработает
    )

    dispatcher.add_handler(conv_handler, 1)
    dispatcher.add_handler(conv_handler2, 2)

    updater.start_polling()

    updater.idle()


if __name__ == '__main__':
    main()
