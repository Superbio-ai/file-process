import json
from typing import List
from io import BytesIO

import pandas as pd
from numpy import number, nan
from pandas.errors import ParserError

from file_process.base import FileProcessorBase
from file_process.constants import PREVIEW_ROWS_COUNT
from file_process.exceptions import NotAllTargetsError, ModelFileValidationVariablesError, DelimiterError, \
    NotSomeTargetsError


class CSVFileProcessor(FileProcessorBase):

    def __init__(self, file: BytesIO, **kwargs):
        delimiter = kwargs.get('delimiter')
        if not delimiter:
            reader = pd.read_csv(file, sep=None, iterator=True, nrows=10)
            delimiter = reader._engine.data.dialect.delimiter  # pylint: disable=protected-access
            file.seek(0)
        self.delimiter = delimiter
        try:
            self.data_df = pd.read_csv(file, sep=self.delimiter)
        except ParserError as exc:
            raise DelimiterError() from exc

    def get_observations(self, rows_number: int = None):
        rows_number = min(10, self.data_df.shape[0]) if rows_number else self.data_df.shape[0]
        return self.data_df.head(rows_number)

    def get_var_names(self):
        return list(self.data_df.columns)

    def get_preview(self):
        var_names = self.get_var_names()
        obs_preview = self.get_observations(PREVIEW_ROWS_COUNT)
        return var_names, self.create_tabular_response(obs_preview), None

    def validate(self):
        pass

    def model_file_validation(self, model_metadata_file: BytesIO):
        reader = json.load(model_metadata_file)
        var_names = set(reader['columns'])
        target_names = set(reader['targets'])
        metadata = reader.get('metadata', {})
        dataset_vars = set(self.data_df.columns)

        all_targets = metadata.get('require_all_targets', 'all')
        if all_targets == 'all':
            difference = target_names - dataset_vars
            if difference:
                raise NotAllTargetsError(difference)
        elif all_targets == 'some':
            are_targets_valid = not target_names or any(elem in dataset_vars for elem in target_names)
            if not are_targets_valid:
                raise NotSomeTargetsError(target_names)
        dataset_diff = dataset_vars - target_names
        var_names_diff = var_names - target_names
        difference = var_names_diff - dataset_diff
        if difference:
            raise ModelFileValidationVariablesError(difference)

    @staticmethod
    def create_tabular_response(data_df: pd.DataFrame) -> List[dict]:
        if data_df is None:
            return []
        numeric_columns = data_df.select_dtypes(include=number).columns
        rows = data_df.replace({nan: None})
        rows[numeric_columns] = rows[numeric_columns].round(2)
        rows = rows.to_dict(orient='records')
        return rows
