import os
import string
import random
import mimetypes
from telegram import Document, File


def gen_id(size, chars=string.ascii_lowercase + string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for x in range(size))


class DocumentBase(Document):

    def __init__(self, dir: str, file_name: str):

        self.local_path = os.path.join(dir, file_name)

        file_id = gen_id(72)
        file_unique_id = gen_id(15)
        file_size = os.path.getsize(self.local_path)
        mime_type = mimetypes.types_map[f".{file_name.split('.')[1]}"]

        super().__init__(file_id=file_id,
                         file_unique_id=file_unique_id,
                         file_name=file_name,
                         mime_type=mime_type,
                         file_size=file_size
                         )


class FileBase(File):

    def __init__(self, document: DocumentBase, file_path: str):

        super().__init__(file_id=document.file_id,
                         file_unique_id=document.file_unique_id,
                         file_name=document.file_name,
                         file_size=document.file_size,
                         file_path=file_path
                         )


class Storage:

    def __init__(self):
        self.storage = {}

    def add(self, file_id, document):
        self.storage[file_id] = document

    def get(self, file_id: str):

        if file_id in self.storage:
            return self.storage[file_id]


storage = Storage()