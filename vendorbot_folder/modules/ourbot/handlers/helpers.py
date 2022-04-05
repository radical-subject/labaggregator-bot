from telegram.ext import CallbackContext


def remove_job_if_exists(name: str, context: CallbackContext) -> bool:
    """
    Удаляет задачу из очереди задач Бота по имени

    :param name: имя задачи. должно быть уникальным.
    :param context:
    :return: False - нет такой задачи, True - задача удалена
    """
    current_jobs = context.job_queue.get_jobs_by_name(name)
    if not current_jobs:
        return False
    for job in current_jobs:
        job.schedule_removal()
    return True


LIST_OF_ADMINS = [336091411, 122267418, 588927967, 47390523]  # tg id's 336091411, 122267418, 588927967, 47390523


def is_admin_chat(chat_id):
    return chat_id in LIST_OF_ADMINS
