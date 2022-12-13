import pytest
from werkzeug.datastructures import FileStorage
import pandas as pd

from file_process.csv import CSVFileProcessor
from file_process.exceptions import ModelFileValidationError, DelimiterError
from tests.test_file_process import get_remote_file_obj, CSV_INPUT_FILES_PATH


class TestCSVFileProcessor:
    original_data_path = f'{CSV_INPUT_FILES_PATH}/original_data.csv'
    MOCK_CONFIGS_PATH = f'{CSV_INPUT_FILES_PATH}/mock_configs'
    MOCK_DATA_PATH = f'{CSV_INPUT_FILES_PATH}/mock_data'

    def test_read_remote_file(self):
        file_obj = get_remote_file_obj(self.original_data_path)
        res = CSVFileProcessor(file_obj)
        assert isinstance(res.data, pd.DataFrame)

    def test_process(self):
        file_obj = get_remote_file_obj(self.original_data_path)
        file_processor = CSVFileProcessor(file_obj)
        var_names, var_preview, obs_preview = file_processor.get_preview()
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
        assert var_preview is None
        assert var_names == ["sepal_length", "sepal_width", "petal_length", "petal_width", "species"]

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
        metadata_file_obj = get_remote_file_obj(data_metadata_path)
        test_file_obj = get_remote_file_obj(new_data_path)
        CSVFileProcessor(test_file_obj).model_file_validation(metadata_file_obj)

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
        metadata_file, metadata_file_obj = get_remote_file_obj(data_metadata_path)
        test_file_obj = get_remote_file_obj(new_data_path)
        with pytest.raises(ModelFileValidationError):
            CSVFileProcessor(test_file_obj).model_file_validation(metadata_file_obj)

    def test_read_file_wrong_delimiter(self):
        test_file_obj = get_remote_file_obj(f'{CSV_INPUT_FILES_PATH}/csv_example.csv')
        with pytest.raises(DelimiterError):
            _ = CSVFileProcessor(test_file_obj, delimiter='.')
