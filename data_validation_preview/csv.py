import json
from typing import List
from io import BytesIO

import pandas as pd
from werkzeug.datastructures import FileStorage

from data_validation_preview.base import FileProcessorBase
from data_validation_preview.exceptions import ModelFileValidationTargetsError, ModelFileValidationVariablesError


class CSVFileProcessor(FileProcessorBase):
    def read_file(self, file, **kwargs):
        if isinstance(file, FileStorage):  # TODO try to get rid of it
            file = file.read()
        file = BytesIO(file)
        read_rows_count = kwargs.get('read_rows_count', 10)
        delimiter = kwargs.get('delimiter', None)
        data = self.read_csv_with_delimiter(file, read_rows_count, delimiter)
        return data

    def read_csv_with_delimiter(self, data_stream, read_rows_count: int, delimiter: str = None):
        if not delimiter:
            reader = pd.read_csv(data_stream, sep=None, iterator=True, nrows=read_rows_count)
            delimiter = reader._engine.data.dialect.delimiter  # pylint: disable=protected-access
            data_stream.seek(0)
        df = pd.read_csv(data_stream, sep=delimiter)
        return df

    def get_preview_data(self, df: pd.DataFrame):
        var_names = list(df.columns)
        obs_preview = df.head(min(10, df.shape[0]))
        return var_names, obs_preview

    def model_file_validation(self, df: pd.DataFrame, model_csv_file: BytesIO, need_target: bool = True):
        # reader = json.load(BytesIO(model_csv_file))
        reader = json.load(model_csv_file)
        var_names = set(reader['columns'])
        target_names = set(reader['targets'])
        metadata = reader.get('metadata', {})
        dataset_vars = set(df.columns)

        if need_target:
            all_targets = metadata.get('require_all_targets', True)
            if all_targets:
                are_targets_valid = not target_names or all(elem in dataset_vars for elem in target_names)
            else:
                are_targets_valid = not target_names or any(elem in dataset_vars for elem in target_names)
            if not are_targets_valid:
                raise ModelFileValidationTargetsError
        are_variables_valid = all(elem in dataset_vars.difference(target_names)
                                  for elem in var_names.difference(target_names))
        if not are_variables_valid:
            raise ModelFileValidationVariablesError

    def process(self, file, model_csv_file: BytesIO = None, **kwargs) \
            -> (None, List[str], dict):
        df = self.read_file(file, **kwargs)
        if model_csv_file:
            self.model_file_validation(df, model_csv_file)
        var_names, obs_preview = self.get_preview_data(df)
        return None, var_names, obs_preview
