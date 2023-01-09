import pytest

from file_process.csv import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.file_processor_factory import FileProcessFactory
from file_process.h5 import H5ADFileProcessor
from tests.test_file_process import get_remote_file_obj
from tests.test_file_process.test_csv import TestCSVFileProcessor
from tests.test_file_process.test_h5 import TestH5ADFileProcessor


class TestFileProcessFactory:
    def test_get_h5_processor(self):
        file_bytes_io = get_remote_file_obj(TestH5ADFileProcessor.path)
        res = FileProcessFactory.get(TestH5ADFileProcessor.path, file_bytes_io)
        assert isinstance(res, H5ADFileProcessor)

    def test_get_h5_processor_wrong_extension(self):
        file_bytes_io = get_remote_file_obj(TestH5ADFileProcessor.path)
        with pytest.raises(WrongExtension):
            _ = FileProcessFactory.get('heart_atlas.h5', file_bytes_io)

    def test_get_csv_processor(self):
        file_bytes_io = get_remote_file_obj(TestCSVFileProcessor.original_data_path)
        res = FileProcessFactory.get(TestCSVFileProcessor.original_data_path, file_bytes_io)
        assert isinstance(res, CSVFileProcessor)
