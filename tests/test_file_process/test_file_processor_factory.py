import pytest

from file_process.csv import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.file_processor_factory import FileProcessFactory
from file_process.h5 import H5ADFileProcessor
from tests.test_file_process import get_remote_file_obj, H5AD_INPUT_FILES_PATH, CSV_INPUT_FILES_PATH


class TestFileProcessFactory:
    def test_get_h5ad_processor(self):
        file = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/heart_sample.h5ad')
        res = FileProcessFactory.get('heart_atlas.h5ad', file)
        assert isinstance(res, H5ADFileProcessor)

    def test_get_csv_processor(self):
        file = get_remote_file_obj(f'{CSV_INPUT_FILES_PATH}/csv_example.csv')
        res = FileProcessFactory.get('heart_atlas.csv', file)
        assert isinstance(res, CSVFileProcessor)

    def test_get_not_implemented_processor(self):
        file = get_remote_file_obj(f'{CSV_INPUT_FILES_PATH}/mock_configs/valid_train_supervised_data_config_1.json')
        with pytest.raises(WrongExtension):
            _ = FileProcessFactory.get('example.json', file)
