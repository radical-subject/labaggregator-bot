import os
import logging
from typing import Dict
from telegram import User, Chat
from .storage import DocumentBase, storage

logger = logging.getLogger(__file__)


BOT_ID = 5000000000

BOT_TOKEN = f'5000000000:BBFJVn-zqLnqQGv_Vrg75aJ5rqppy410rm0'

CHAT_ID = 1


class UserBase(User):

    def __init__(self, id: int = CHAT_ID):
        logging.info(f'create UserBase={id}')
        super().__init__(
            id,  # id
            'FirstName',  # first_name
            False,  # is_bot
            'LastName',  # last_name
            f'user{id}',  # username
            'ru'  # language_code
        )


class ChatBase(Chat):

    def __init__(self, user: UserBase):
        super().__init__(
            user.id,
            'private',  # type
            None,  # title
            user.username,
            user.first_name,
            user.last_name,
        )


virtual_bot = User(
    BOT_ID,
    'PythonTelegramUnitTestBot',  # first_name
    True,  # is_bot
    None,  # last_name
    'PTUTestBot',  # username
    None,  # language_code
    True,  # can_join_groups
    False,  # can_read_all_group_messages
    False  # supports_inline_queries
)


class Tester:

    def __init__(self, core, user, chat):
        self.core = core
        self.user = user
        self.chat = chat

    def init_dialog(self):
        self.send_command('/start')
        self.clean_dialog()

    def clean_dialog(self):
        while self.get_message(timeout=0.5):
            pass

    def send_message(self, text: str) -> None:

        self.core.user_send(virtual_bot.id,
                            user_from=self.user.to_dict(),
                            chat=self.chat.to_dict(),
                            text=text
                            )

    def send_command(self, command: str) -> None:

        self.core.user_send_command(virtual_bot.id,
                                    user_from=self.user.to_dict(),
                                    chat=self.chat.to_dict(),
                                    command=command
                                    )

    def get_message(self, timeout=2.0) -> Dict:
        messages = self.core.get_updates(self.user.id, timeout)
        if messages:
            return messages[0]['message']

    def send_file(self, dir: str, file_name: str) -> None:

        try:
            document = DocumentBase(dir, file_name)
            storage.add(document['file_id'], document)

            self.core.user_send(virtual_bot.id,
                                user_from=self.user.to_dict(),
                                chat=self.chat.to_dict(),
                                document=document.to_dict()
                                )

        except Exception as err:
            logger.error(err)

