import os
from modules.bot.bot import BotObject

from setup_logger import setup_logger


def main():
    token = os.getenv('BOT_TOKEN')
    if not token:
        raise Exception('set BOT_TOKEN variable')

    setup_logger()
    bot = BotObject(token)
    bot.start()


if __name__ == '__main__':
    main()
