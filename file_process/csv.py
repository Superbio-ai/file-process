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
        try:
            df = pd.read_csv(data_stream, sep=delimiter)
        except ParserError as exc:
            raise DelimiterError() from exc
        return df

    def get_preview_data(self, df: pd.DataFrame):
        var_names = list(df.columns)
        obs_preview = df.head(min(10, df.shape[0]))
        return var_names, obs_preview

    def model_file_validation(self, df: pd.DataFrame, model_metadata_file: BytesIO, need_target: bool = True):
        reader = json.load(model_metadata_file)
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

    def process(self, file, model_metadata_file: BytesIO = None, **kwargs) -> (List[str], None, pd.DataFrame):
        df = self.read_file(file, **kwargs)
        if model_metadata_file:
            self.model_file_validation(df, model_metadata_file)
        var_names, obs_preview = self.get_preview_data(df)
        return var_names, None, obs_preview

    def create_tabular_response(self, data_df: pd.DataFrame) -> List[dict]:
        if data_df is None:
            return []
        numeric_columns = data_df.select_dtypes(include=number).columns
        rows = data_df.replace({nan: None})
        rows[numeric_columns] = rows[numeric_columns].round(2)
        rows = rows.to_dict(orient='records')
        return rows
