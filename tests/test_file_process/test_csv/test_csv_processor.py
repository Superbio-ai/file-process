from io import BytesIO

import pytest
from file_process.exceptions import DelimiterError
from numpy import nan
import pandas as pd

from file_process.csv.csv_processor import CSVFileProcessor
from tests.test_file_process import CSV_INPUT_FILES_PATH, get_remote_file_obj


class TestCSVFileProcessor:
    original_data_path = f'{CSV_INPUT_FILES_PATH}/original_data.csv'
    single_column_data_path = f'{CSV_INPUT_FILES_PATH}/single_column_data.csv'
    tsv_data_path = f'{CSV_INPUT_FILES_PATH}/tsv_example.tsv'
    MOCK_CONFIGS_PATH = f'{CSV_INPUT_FILES_PATH}/mock_configs'
    MOCK_DATA_PATH = f'{CSV_INPUT_FILES_PATH}/mock_data'

    def _get_file_and_remote_file_obj(self, path: str):
        file = open(path, 'rb')
        file_obj = BytesIO(file.read())
        return file, file_obj

    def test_read_file(self):
        file_bytes_io = get_remote_file_obj(self.original_data_path)
        res = CSVFileProcessor(file_bytes_io)
        assert isinstance(res.data_df, pd.DataFrame)

    def test_get_preview(self):
        file_bytes_io = get_remote_file_obj(self.original_data_path)
        file_processor = CSVFileProcessor(file_bytes_io)
        var_names, obs_preview, var_preview, _ = file_processor.get_preview()
        assert obs_preview == [
            {"sepal_length": 5.1, "sepal_width": 3.5, "petal_length": 1.4, "petal_width": 0.2, "species": "setosa"},
            {"sepal_length": 4.9, "sepal_width": 3.0, "petal_length": 1.4, "petal_width": 0.2, "species": "setosa"},
            {"sepal_length": 4.7, "sepal_width": 3.2, "petal_length": 1.3, "petal_width": 0.2, "species": "setosa"},
            {"sepal_length": 4.6, "sepal_width": 3.1, "petal_length": 1.5, "petal_width": 0.2, "species": "setosa"},
            {"sepal_length": 5.0, "sepal_width": 3.6, "petal_length": 1.4, "petal_width": 0.2, "species": "setosa"},
            {"sepal_length": 5.4, "sepal_width": 3.9, "petal_length": 1.7, "petal_width": 0.4, "species": "setosa"},
            {"sepal_length": 4.6, "sepal_width": 3.4, "petal_length": 1.4, "petal_width": 0.3, "species": "setosa"},
            {"sepal_length": 5.0, "sepal_width": 3.4, "petal_length": 1.5, "petal_width": 0.2, "species": "setosa"},
            {"sepal_length": 4.4, "sepal_width": 2.9, "petal_length": 1.4, "petal_width": 0.2, "species": "setosa"},
            {"sepal_length": 4.9, "sepal_width": 3.1, "petal_length": 1.5, "petal_width": 0.1, "species": "setosa"}
        ]
        assert var_preview == None
        assert var_names == ["sepal_length", "sepal_width", "petal_length", "petal_width", "species"]

    def test_get_preview_single_column(self):
        file_bytes_io = get_remote_file_obj(self.single_column_data_path)
        file_processor = CSVFileProcessor(file_bytes_io)
        var_names, obs_preview, var_preview, _ = file_processor.get_preview()
        assert obs_preview == [
            {"sepal_length": 5.1},
            {"sepal_length": 4.9},
            {"sepal_length": 4.7},
            {"sepal_length": 4.6},
            {"sepal_length": 5.0},
            {"sepal_length": 5.4},
            {"sepal_length": 4.6},
        ]
        assert var_preview == None
        assert var_names == ["sepal_length"]

    def test_get_preview_file_with_nans(self):
        file_bytes_io = get_remote_file_obj(f'{CSV_INPUT_FILES_PATH}/follicular_obs_sample.csv')
        var_names, obs_preview, var_preview, _ = CSVFileProcessor(file_bytes_io).get_preview()
        for item in obs_preview:
            for value in item.values():
                assert value is not nan

    def test_read_file_wrong_delimiter(self):
        file_bytes_io = get_remote_file_obj(f'{CSV_INPUT_FILES_PATH}/csv_example.csv')
        with pytest.raises(DelimiterError):
            _ = CSVFileProcessor(file_bytes_io, delimiter='.')

    def test_tsv_preview(self):
        file_bytes_io = get_remote_file_obj(self.tsv_data_path)
        file_processor = CSVFileProcessor(file_bytes_io)
        var_names, obs_preview, var_preview, _ = file_processor.get_preview()
        assert obs_preview == [
            {"FID": 0, "IID": 0, "cov1": 0.09, "cov2": 0.69, "target_label": 1.2},
            {"FID": 1, "IID": 1, "cov1": 0.24, "cov2": 0.54, "target_label": 1.1},
            {"FID": 2, "IID": 2, "cov1": 0.31, "cov2": 0.01, "target_label": 0.3},
        ]
        assert var_preview == None
        assert var_names == ["FID", "IID", "cov1", "cov2", "target_label"]

