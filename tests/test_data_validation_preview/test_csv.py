import pytest
from werkzeug.datastructures import FileStorage
import pandas as pd

from common import create_tabular_csv_response
from data_validation_preview.csv import CSVFileProcessor
from data_validation_preview.exceptions import ModelFileValidationError
from tests.test_data_validation_preview import INPUT_FILES_PATH


class TestCSVFileProcessor:
    MAIN_PATH = f'{INPUT_FILES_PATH}/tabular_csv'
    original_data_path = f'{MAIN_PATH}/original_data.csv'
    MOCK_CONFIGS_PATH = f'{MAIN_PATH}/mock_configs'
    MOCK_DATA_PATH = f'{MAIN_PATH}/mock_data'

    def _get_file_and_remote_file_obj(self, path: str):
        with open(path, 'rb') as file:
            file_obj = file.read()
            return file, file_obj

    def test_read_local_file(self):
        with open(self.original_data_path, 'rb') as file:
            file_obj = FileStorage(file)
            res = CSVFileProcessor().read_file(file_obj)
            assert isinstance(res, pd.DataFrame)

    def test_read_remote_file(self):
        file, file_obj = self._get_file_and_remote_file_obj(self.original_data_path)
        res = CSVFileProcessor().read_file(file_obj)
        assert isinstance(res, pd.DataFrame)
        file.close()

    def test_process(self):
        file, file_obj = self._get_file_and_remote_file_obj(self.original_data_path)
        var_names, var_preview, obs_preview = CSVFileProcessor().process(file_obj)
        obs_preview = create_tabular_csv_response(obs_preview)
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
        assert var_preview == []
        original_var_names = ["sepal_length", "sepal_width", "petal_length", "petal_width", "species"]
        assert var_names == ["sepal_length", "sepal_width", "petal_length", "petal_width", "species"]
        assert all(var in var_names for var in original_var_names)
        file.close()

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

    @pytest.mark.parametrize('config_csv, data_csv', valid_tuples)
    def tst_model_file_validation_with_csv(self, config_csv, data_csv):
        new_data_path = f'{self.MOCK_DATA_PATH}/{data_csv}'
        data_metadata_path = f'{self.MOCK_CONFIGS_PATH}/{config_csv}'
        metadata_file, metadata_file_obj = self._get_file_and_remote_file_obj(data_metadata_path)
        test_file, test_file_obj = self._get_file_and_remote_file_obj(new_data_path)
        df = CSVFileProcessor().read_file(test_file_obj)
        is_valid = CSVFileProcessor().model_file_validation(df, metadata_file_obj)
        metadata_file.close()
        test_file.close()

    invalid_tuples = [
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_1.json'),
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_2.json'),
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_3.json'),
        ('sometimes_valid_new_data.csv', 'valid_train_supervised_data_config_5.json'),
        ('sometimes_valid_new_data.csv', 'valid_usupervised_data_config_1.json'),
        ('sometimes_valid_new_data.csv', 'valid_usupervised_data_config_2.json'),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_1.json'),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_2.json'),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_3.json'),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_4.json'),
        ('invalid_new_data.csv', 'valid_train_supervised_data_config_5.json'),
        ('invalid_new_data.csv', 'valid_usupervised_data_config_1.json'),
        ('invalid_new_data.csv', 'valid_usupervised_data_config_2.json'),
    ]

    @pytest.mark.parametrize('config_csv, data_csv', invalid_tuples)
    def tst_model_file_validation_with_csv_error(self, config_csv, data_csv):
        new_data_path = f'{self.MOCK_DATA_PATH}/{data_csv}'
        data_metadata_path = f'{self.MOCK_CONFIGS_PATH}/{config_csv}'
        metadata_file, metadata_file_obj = self._get_file_and_remote_file_obj(data_metadata_path)
        test_file, test_file_obj = self._get_file_and_remote_file_obj(new_data_path)
        df = CSVFileProcessor().read_file(test_file_obj)
        with pytest.raises(ModelFileValidationError):
            is_valid = CSVFileProcessor().model_file_validation(df, metadata_file_obj)
        metadata_file.close()
        test_file.close()
