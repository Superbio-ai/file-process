from io import BytesIO

import pytest

from file_process.csv.csv_processor import CSVFileProcessor
from file_process.csv.csv_validator import CSVValidator
from file_process.exceptions import NotAllTargetsError, ModelFileValidationVariablesError, CustomValidationException
from tests.test_file_process import CSV_INPUT_FILES_PATH, get_remote_file_obj


class TestCSVValidator:
    original_data_path = f'{CSV_INPUT_FILES_PATH}/original_data.csv'
    MOCK_CONFIGS_PATH = f'{CSV_INPUT_FILES_PATH}/mock_configs'
    MOCK_DATA_PATH = f'{CSV_INPUT_FILES_PATH}/mock_data'

    def _get_file_and_remote_file_obj(self, path: str):
        file = open(path, 'rb')
        file_obj = BytesIO(file.read())
        return file, file_obj

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
        _ = CSVFileProcessor(file_bytes_io).validate(metadata_file_bytes_io)

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
            _ = CSVFileProcessor(file_bytes_io).validate(metadata_file_bytes_io)

    validation_dicts_for_valid_tests = [
        {
            "columnsList": [{"name": "sepal_length"}, {"name": "sepal_width"}, {"name": "petal_length"},
                            {"name": "petal_width"}, {"name": "species"}]
        }, {
            "columnsList": [{"name": "sepal_width"}, {"name": "sepal_length"}, {"name": "petal_length"},
                            {"name": "petal_width"}, {"name": "species"}],
            "preserveOrder": False
        }, {
            "columnsList": [{"name": "sepal_length"}, {"name": "sepal_width"}, {"name": "petal_length"}],
            "allowOtherColumns": True
        }
    ]

    @pytest.mark.parametrize('validation_dict', validation_dicts_for_valid_tests)
    def test_validate_columnns_valid(self, validation_dict):
        file_bytes_io = get_remote_file_obj(self.original_data_path)
        processor = CSVFileProcessor(file_bytes_io)

        validation_rules = {"columns": validation_dict}
        validator = CSVValidator(processor.data_df, validation_rules)
        validator._validate_column_names()

    validation_dicts_for_invalid_tests = [
        {
            "columnsList": [{"name": "some_name"}, {"name": "sepal_width"}, {"name": "petal_length"},
                            {"name": "petal_width"}, {"name": "species"}]
        }, {
            "columnsList": [{"name": "sepal_width"}, {"name": "sepal_length"}, {"name": "petal_length"},
                            {"name": "petal_width"}, {"name": "species"}],
            "preserveOrder": True
        }, {
            "columnsList": [{"name": "sepal_length"}, {"name": "sepal_width"}, {"name": "petal_length"}],
            "allowOtherColumns": False
        }
    ]

    @pytest.mark.parametrize('validation_dict', validation_dicts_for_invalid_tests)
    def test_validate_columnns_invalid(self, validation_dict):
        file_bytes_io = get_remote_file_obj(self.original_data_path)
        processor = CSVFileProcessor(file_bytes_io)

        validation_rules = {"columns": validation_dict}
        validator = CSVValidator(processor.data_df, validation_rules)
        with pytest.raises(CustomValidationException):
            validator._validate_column_names()
