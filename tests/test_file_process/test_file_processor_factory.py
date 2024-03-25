import pytest

from file_process.csv.csv_processor import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.file_processor_factory import FileProcessFactory
from file_process.h5ad.h5ad_processor import H5ADFileProcessor
from file_process.txt.txt_processor import TxtFileProcessor
from tests.test_file_process import get_remote_file_obj
from tests.test_file_process.test_csv.test_csv_processor import TestCSVFileProcessor
from tests.test_file_process.test_h5ad.test_h5ad_processor import TestH5ADFileProcessor
from tests.test_file_process.test_txt.test_txt_processor import TestTxtFileProcessor


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

    def test_get_csv_processor_tsv_file(self):
        file_bytes_io = get_remote_file_obj(TestCSVFileProcessor.tsv_data_path)
        res = FileProcessFactory.get(TestCSVFileProcessor.tsv_data_path, file_bytes_io)
        assert isinstance(res, CSVFileProcessor)

    def test_get_txt_processor(self):
        file_bytes_io = get_remote_file_obj(TestTxtFileProcessor.original_data_path)
        res = FileProcessFactory.get(TestTxtFileProcessor.original_data_path, file_bytes_io)
        assert isinstance(res, TxtFileProcessor)

    def test_get_txt_processor_pdb_file(self):
        file_bytes_io = get_remote_file_obj(TestTxtFileProcessor.pdb_data_path)
        res = FileProcessFactory.get(TestTxtFileProcessor.pdb_data_path, file_bytes_io)
        assert isinstance(res, TxtFileProcessor)

    def test_get_txt_processor_fasta_file(self):
        file_bytes_io = get_remote_file_obj(TestTxtFileProcessor.fasta_data_path)
        res = FileProcessFactory.get(TestTxtFileProcessor.fasta_data_path, file_bytes_io)
        assert isinstance(res, TxtFileProcessor)

    def test_validate_wrong_extension(self):
        filename = 'fastq.fastq'
        with pytest.raises(WrongExtension):
            _ = FileProcessFactory.validate_extension(filename)

    def test_validate_extension(self):
        filename = 'fastq.h5ad'
        FileProcessFactory.validate_extension(filename)
