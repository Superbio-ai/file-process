from typing import List, Optional
from io import BytesIO
import csv
import _csv

import pandas as pd
from numpy import number, nan
from pandas.errors import ParserError

from file_process.base import FileProcessorBase
from file_process.constants import PREVIEW_ROWS_COUNT
from file_process.csv.csv_validator import CSVValidator
from file_process.csv.schemas import SbioModelDataForCsv
from file_process.exceptions import DelimiterError


class CSVFileProcessor(FileProcessorBase):

    def __init__(self, file: BytesIO, **kwargs):
        delimiter = kwargs.get('delimiter')
        if not delimiter:
            sniffer = csv.Sniffer()
            data = file.read(4096)
            try:
                delimiter = sniffer.sniff(str(data, encoding='utf-8')).delimiter
            except _csv.Error:
                delimiter = ","
            file.seek(0)
        self.delimiter = delimiter
        try:
            self.data_df = pd.read_csv(file, sep=self.delimiter)
        except ParserError as exc:
            raise DelimiterError() from exc

    def validate(self, model_metadata_file: Optional[BytesIO] = None, validation_rules: Optional[dict] = None):
        model_data = SbioModelDataForCsv(model_metadata_file) if model_metadata_file else None
        validator = CSVValidator(self.data_df, validation_rules, model_data)
        validator()

    def get_observations(self, rows_number: int = None):
        rows_number = min(10, self.data_df.shape[0]) if rows_number else self.data_df.shape[0]
        return self.data_df.head(rows_number)

    def get_var_names(self):
        return list(self.data_df.columns)

    def get_preview(self):
        var_names = self.get_var_names()
        obs_preview = self.get_observations(PREVIEW_ROWS_COUNT)
        return var_names, None, self.create_tabular_response(obs_preview), None, None

    @staticmethod
    def create_tabular_response(data_df: pd.DataFrame) -> List[dict]:
        if data_df is None:
            return []
        numeric_columns = data_df.select_dtypes(include=number).columns
        rows = data_df.replace({nan: None})
        rows[numeric_columns] = rows[numeric_columns].round(2)
        rows = rows.to_dict(orient='records')
        return rows
