from file_process.csv import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.h5 import H5FileProcessor


class FileProcessFactory:  # pylint: disable=too-few-public-methods
    @classmethod
    def get(cls, filename: str):
        if filename.endswith('.h5') or filename.endswith('.h5ad'):
            return H5FileProcessor()
        if filename.endswith('.csv'):
            return CSVFileProcessor()
        raise WrongExtension
