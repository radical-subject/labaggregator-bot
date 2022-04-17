import traceback
import logging
logger = logging.getLogger(__name__)

from functools import wraps
from modules.ourbot.handlers.helpers import is_admin_chat


def log_errors(f):
    """
    #TODO кажется в Bot есть обработчик ошибок.
    # Также есть https://github.com/python-telegram-bot/python-telegram-bot/tests/test_dispatcher.py#L153
    # add_error_handler(self.error_handler)
    декоратор @log_errors который импортируется в файл и ставится перед функциями,
    которым нужен le-дебаг
    """
    def inner(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as e:
            logger.error(traceback.format_exc())
            raise e
    return inner


## ADMINISTRATION decorator
def is_admin(func):
    """
    декоратор запрещающий пользование функциями перед которыми он стоит 
    аккаунтами с id которые не совпадают с теми, что в списке выше
    """
    @wraps(func)
    def wrapped(self, update, context, *args, **kwargs):
        user_id = update.effective_user.id
        if not is_admin_chat(user_id):
            logger.error(f"Unauthorized access denied for {user_id}.")
            return
        return func(self, update, context, *args, **kwargs)
    return wrapped
