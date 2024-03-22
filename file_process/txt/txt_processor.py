from typing import List, Optional
from io import BytesIO
import csv

import pandas as pd
from numpy import number, nan
from pandas.errors import ParserError

from file_process.base import FileProcessorBase
from file_process.constants import PREVIEW_ROWS_COUNT
from file_process.csv.csv_validator import CSVValidator
from file_process.csv.schemas import SbioModelDataForCsv
from file_process.exceptions import DelimiterError


class TxtFileProcessor(FileProcessorBase):

    def __init__(self, file: BytesIO, **kwargs):
        self.txt_data = file.read()

    def get_preview(self):
        return None, None, None, self.txt_data

    def validate(self, model_metadata_file: Optional[BytesIO] = None, validation_rules: Optional[dict] = None):
        pass
