import json
from typing import List
from io import BytesIO

import pandas as pd
from numpy import number, nan
from pandas.errors import ParserError
from werkzeug.datastructures import FileStorage

from file_process.base import FileProcessorBase
from file_process.exceptions import ModelFileValidationTargetsError, ModelFileValidationVariablesError, DelimiterError


class CSVFileProcessor(FileProcessorBase):

    def __init__(self, file, **kwargs):
        if isinstance(file, FileStorage):  # TODO try to get rid of it
            file = file.read()
            file = BytesIO(file)
        read_rows_count = kwargs.get('read_rows_count', 10)
        delimiter = kwargs.get('delimiter', None)
        self.data = self._read_csv_with_delimiter(file, read_rows_count, delimiter)

    @staticmethod
    def _read_csv_with_delimiter(data_stream, read_rows_count: int, delimiter: str = None):
        if not delimiter:
            reader = pd.read_csv(data_stream, sep=None, iterator=True, nrows=read_rows_count)
            delimiter = reader._engine.data.dialect.delimiter  # pylint: disable=protected-access
            data_stream.seek(0)
        try:
            df = pd.read_csv(data_stream, sep=delimiter)
        except ParserError as exc:
            raise DelimiterError() from exc
        return df

    def get_obs(self):
        return self.data.head(min(10, self.data.shape[0]))

    def get_var_names(self):
        return list(self.data.columns)

    def get_preview(self):
        var_names = self.get_var_names()

        obs = self.get_obs()
        obs_preview = self._create_tabular_response(obs)

        return var_names, obs_preview, None

    def model_file_validation(self, model_metadata_file: BytesIO, need_target: bool = True):
        reader = json.load(model_metadata_file)
        var_names = set(reader['columns'])
        target_names = set(reader['targets'])
        metadata = reader.get('metadata', {})
        dataset_vars = set(self.data.columns)

        are_variables_valid = all(elem in dataset_vars.difference(target_names)
                                  for elem in var_names.difference(target_names))
        if not are_variables_valid:
            raise ModelFileValidationVariablesError

        if not need_target:
            return
        all_targets = metadata.get('require_all_targets', True)
        if all_targets:
            are_targets_valid = not target_names or all(elem in dataset_vars for elem in target_names)
            if not are_targets_valid:
                raise ModelFileValidationTargetsError
        else:
            are_targets_valid = not target_names or any(elem in dataset_vars for elem in target_names)
            if not are_targets_valid:
                raise ModelFileValidationTargetsError

    def validate(self):
        pass

    @staticmethod
    def _create_tabular_response(data_df: pd.DataFrame) -> List[dict]:
        if data_df is None:
            return []
        numeric_columns = data_df.select_dtypes(include=number).columns
        rows = data_df.replace({nan: None})
        rows[numeric_columns] = rows[numeric_columns].round(2)
        rows = rows.to_dict(orient='records')
        return rows
