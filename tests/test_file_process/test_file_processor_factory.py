from io import BytesIO

import pytest

from file_process.csv import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.file_processor_factory import FileProcessFactory
from file_process.h5 import H5FileProcessor
from tests.test_file_process import INPUT_FILES_PATH


class TestFileProcessFactory:
    def _get_file_and_remote_file_obj(self, path: str):
        file = open(path, 'rb')
        file_obj = BytesIO(file.read())
        return file, file_obj

    def test_get_h5_processor(self):
        file, file_obj = self._get_file_and_remote_file_obj(f'{INPUT_FILES_PATH}/heart_sample.h5ad')
        res = FileProcessFactory.get('heart_atlas.h5ad', file_obj)
        assert isinstance(res, H5FileProcessor)
        file.close()

    def test_get_h5ad_processor(self):
        res = FileProcessFactory.get('heart_atlas.h5')
        assert isinstance(res, H5FileProcessor)

    def test_get_csv_processor(self):
        res = FileProcessFactory.get('heart_atlas.csv')
        assert isinstance(res, CSVFileProcessor)

    def test_get_not_implemented_processor(self):
        with pytest.raises(WrongExtension):
            _ = FileProcessFactory.get('heart_atlas.doc')
