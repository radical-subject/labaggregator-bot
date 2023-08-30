import sys
import os
from logging import getLogger, getLevelName, StreamHandler, Formatter
from logging.handlers import RotatingFileHandler


# Set up the logging configuration
log_filename = 'bot.log'
log_max_size = 10 * 1024 * 1024  # 10 MB
log_backup_count = 5  # Number of backup log files to keep: bot.log, bot.log.1, ..., bot.log.5

logs_dir = os.path.join(os.path.dirname(__file__), os.path.pardir)  # current file parent folder
log_file_path = os.path.normpath(os.path.join(logs_dir, log_filename))


def setup_logger():
    """
    Можно задать переменной окружения LOG_LEVEL.
    Значения: DEBUG, INFO и др.
    """
    log = getLogger()
    log_level = getLevelName(os.getenv("LOG_LEVEL", "INFO").upper())

    log.setLevel(log_level)

    #if not os.path.isdir(logs_dir):   # не надо т.к. в директории репозитория
    #    os.makedirs(logs_dir)

    # Create a rotating file handler
    fh = RotatingFileHandler(
        log_file_path, maxBytes=log_max_size, backupCount=log_backup_count, mode='w'
    )

    fh.setLevel(log_level)

    formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    fh.setFormatter(formatter)
    log.addHandler(fh)

    console = sys.stderr
    if console is not None:
        console = StreamHandler(console)
        console.setLevel(log_level)
        console.setFormatter(formatter)
        log.addHandler(console)

    log.info(f"log level={log_level}")
