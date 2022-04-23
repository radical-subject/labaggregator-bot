import os
import sys
from logging import getLogger, getLevelName, root, Formatter, StreamHandler

from modules.token_extractor import token
from modules.bot.bot import BotObject

logger = getLogger()


def setup_logger():
    logger.setLevel(getLevelName(os.getenv("LOG_LEVEL", "INFO").upper()))

    formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    console = sys.stderr
    if console is not None:
        console = StreamHandler(console)
        console.setFormatter(formatter)
        logger.addHandler(console)

    logger.info('log level=' + getLevelName(root.level))


def main():

    setup_logger()
    bot = BotObject(token)
    bot.start()


if __name__ == '__main__':
    main()
