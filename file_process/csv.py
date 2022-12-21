import json
from typing import List
from io import BytesIO

import pandas as pd
from numpy import number, nan
from pandas.errors import ParserError

from file_process.base import FileProcessorBase
from file_process.exceptions import NotAllTargetsError, ModelFileValidationVariablesError, DelimiterError, \
    NotSomeTargetsError


class CSVFileProcessor(FileProcessorBase):

    def __init__(self, file: BytesIO, **kwargs):
        self.delimiter = kwargs.get('delimiter') or self._get_delimiter(file)
        try:
            self.data_df = pd.read_csv(file, sep=self.delimiter)
        except ParserError as exc:
            raise DelimiterError() from exc

    def _get_delimiter(self, file: BytesIO):
        reader = pd.read_csv(file, sep=None, iterator=True, nrows=10)
        delimiter = reader._engine.data.dialect.delimiter  # pylint: disable=protected-access
        file.seek(0)
        return delimiter

    def get_obs(self):
        return self.data_df.head(min(10, self.data_df.shape[0]))

    def get_var_names(self):
        return list(self.data_df.columns)

    def get_preview(self):
        var_names = self.get_var_names()
        obs_preview = self.get_obs()
        return var_names, self.create_tabular_response(obs_preview), None

    def validate(self):
        pass

    def model_file_validation(self, model_metadata_file: BytesIO):
        reader = json.load(model_metadata_file)
        var_names = set(reader['columns'])
        target_names = set(reader['targets'])
        metadata = reader.get('metadata', {})
        dataset_vars = set(self.data_df.columns)

        all_targets = metadata.get('require_all_targets', True)
        if all_targets == 'all':
            are_targets_valid = not target_names or all(elem in dataset_vars for elem in target_names)
            if not are_targets_valid:
                raise NotAllTargetsError
        elif all_targets == 'some':
            are_targets_valid = not target_names or any(elem in dataset_vars for elem in target_names)
            if not are_targets_valid:
                raise NotSomeTargetsError

        are_variables_valid = all(elem in dataset_vars.difference(target_names)
                                  for elem in var_names.difference(target_names))
        if not are_variables_valid:
            raise ModelFileValidationVariablesError

    @staticmethod
    def create_tabular_response(data_df: pd.DataFrame) -> List[dict]:
        if data_df is None:
            return []
        numeric_columns = data_df.select_dtypes(include=number).columns
        rows = data_df.replace({nan: None})
        rows[numeric_columns] = rows[numeric_columns].round(2)
        rows = rows.to_dict(orient='records')
        return rows
