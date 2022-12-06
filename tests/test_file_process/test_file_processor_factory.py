import pytest

from file_process.csv import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.file_processor_factory import FileProcessFactory
from file_process.h5 import H5FileProcessor


class TestFileProcessFactory:
    def test_get_h5_processor(self):
        res = FileProcessFactory.get('heart_atlas.h5ad')
        assert isinstance(res, H5FileProcessor)

    def test_get_h5ad_processor(self):
        res = FileProcessFactory.get('heart_atlas.h5')
        assert isinstance(res, H5FileProcessor)

    def test_get_csv_processor(self):
        res = FileProcessFactory.get('heart_atlas.csv')
        assert isinstance(res, CSVFileProcessor)

    def test_get_not_implemented_processor(self):
        with pytest.raises(WrongExtension):
            _ = FileProcessFactory.get('heart_atlas.doc')
