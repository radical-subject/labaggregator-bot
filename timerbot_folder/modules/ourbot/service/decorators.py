import logging, traceback

logger = logging.getLogger(__name__)


def log_errors(f):
    """
    декоратор @log_errors который импортируется в файл и ставится перед функциями,
    которым нужен le-дебаг
    """
    def inner(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as e:
            logger.error(f'{e}')
            e_traceback = traceback.format_exception(e.__class__, e, e.__traceback__)
            traceback_lines = []
            for line in [line.rstrip('\n') for line in e_traceback]:
                traceback_lines.extend(line.splitlines())
            error_message = '\n'.join(traceback_lines)
            logger.error(error_message)
    return inner


## ADMINISTRATION decorator
from functools import wraps
LIST_OF_ADMINS = [336091411, 122267418, 588927967] #tg id's 336091411, 122267418, 588927967
def restricted(func):
    """
    декоратор запрещающий пользование функциями перед которыми он стоит 
    аккаунтами с id которые не совпадают с теми, что в списке выше
    """
    @wraps(func)
    def wrapped(self, update, context, *args, **kwargs):
        user_id = update.effective_user.id
        if user_id not in LIST_OF_ADMINS:
            print("Unauthorized access denied for {}.".format(user_id))
            return
        return func(self, update, context, *args, **kwargs)
    return wrapped
