from io import BytesIO

import pytest
from numpy import nan
import pandas as pd

from file_process.csv.csv_processor import CSVFileProcessor
from file_process.exceptions import DelimiterError, NotAllTargetsError, ModelFileValidationVariablesError
from tests.test_file_process import CSV_INPUT_FILES_PATH, get_remote_file_obj


class TestCSVFileProcessor:
    original_data_path = f'{CSV_INPUT_FILES_PATH}/original_data.csv'
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
        var_names, obs_preview, var_preview = file_processor.get_preview()
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

    def test_get_preview_file_with_nans(self):
        file_bytes_io = get_remote_file_obj(f'{CSV_INPUT_FILES_PATH}/follicular_obs_sample.csv')
        var_names, obs_preview, var_preview = CSVFileProcessor(file_bytes_io).get_preview()
        for item in obs_preview:
            for value in item.values():
                assert value is not nan

    valid_tuples = [
        ('valid_new_data.csv', 'valid_train_supervised_data_config_1.json'),
        ('valid_new_data.csv', 'valid_train_supervised_data_config_2.json'),
        ('valid_new_data.csv', 'valid_train_supervised_data_config_3.json'),
        ('valid_new_data.csv', 'valid_train_supervised_data_config_4.json'),
        ('valid_new_data.csv', 'valid_train_supervised_data_config_5.json'),
        ('valid_new_data.csv', 'valid_usupervised_data_config_1.json'),
        ('valid_new_data.csv', 'valid_usupervised_data_config_2.json'),
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_4.json')
    ]

    @pytest.mark.parametrize('data_csv, config_csv', valid_tuples)
    def test_model_file_validation_with_csv(self, data_csv, config_csv):
        file_bytes_io = get_remote_file_obj(f'{self.MOCK_DATA_PATH}/{data_csv}')
        metadata_file_bytes_io = get_remote_file_obj(f'{self.MOCK_CONFIGS_PATH}/{config_csv}')
        _ = CSVFileProcessor(file_bytes_io).model_file_validation(metadata_file_bytes_io)

    invalid_tuples = [
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_1.json', NotAllTargetsError),
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_2.json', NotAllTargetsError),
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_3.json', NotAllTargetsError),
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_5.json', NotAllTargetsError),
        ('sometimes_valid_new_data.csv', 'valid_usupervised_data_config_1.json', ModelFileValidationVariablesError),
        ('sometimes_valid_new_data.csv', 'valid_usupervised_data_config_2.json', ModelFileValidationVariablesError),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_1.json', ModelFileValidationVariablesError),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_2.json', ModelFileValidationVariablesError),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_3.json', ModelFileValidationVariablesError),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_4.json', ModelFileValidationVariablesError),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_5.json', ModelFileValidationVariablesError),
        ('invalid_new_data.csv', 'valid_usupervised_data_config_1.json', ModelFileValidationVariablesError),
        ('invalid_new_data.csv', 'valid_usupervised_data_config_2.json', ModelFileValidationVariablesError),
    ]

    @pytest.mark.parametrize('data_csv, config_csv, exception', invalid_tuples)
    def test_model_file_validation_with_csv_error(self, data_csv, config_csv, exception):
        file_bytes_io = get_remote_file_obj(f'{self.MOCK_DATA_PATH}/{data_csv}')
        metadata_file_bytes_io = get_remote_file_obj(f'{self.MOCK_CONFIGS_PATH}/{config_csv}')
        with pytest.raises(exception):
            CSVFileProcessor(file_bytes_io).model_file_validation(metadata_file_bytes_io)

    def test_read_file_wrong_delimiter(self):
        file_bytes_io = get_remote_file_obj(f'{CSV_INPUT_FILES_PATH}/csv_example.csv')
        with pytest.raises(DelimiterError):
            _ = CSVFileProcessor(file_bytes_io, delimiter='.')
