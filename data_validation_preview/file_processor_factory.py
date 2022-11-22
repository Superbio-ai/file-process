from data_validation_preview.csv import CSVFileProcessor
from data_validation_preview.exceptions import WrongExtension
from data_validation_preview.h5 import H5FileProcessor


class FileProcessFactory:  # pylint: disable=too-few-public-methods
    @classmethod
    def get_processor(cls, filename: str):
        if filename.endswith('.h5') or filename.endswith('.h5ad'):
            return H5FileProcessor()
        if filename.endswith('.csv'):
            return CSVFileProcessor()
        raise WrongExtension
