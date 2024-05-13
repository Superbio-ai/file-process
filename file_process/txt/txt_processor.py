from typing import Optional
from io import BytesIO

from file_process.base import FileProcessorBase


class TxtFileProcessor(FileProcessorBase):

    def __init__(self, file: BytesIO, **_):
        self.txt_data = str(file.read(), encoding='utf-8')

    def get_preview(self):
        return None, None, None, None, self.txt_data

    def validate(self, model_metadata_file: Optional[BytesIO] = None, validation_rules: Optional[dict] = None):
        pass
