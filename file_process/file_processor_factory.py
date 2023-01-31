from io import BytesIO

from file_process.csv import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.h5 import H5ADFileProcessor


class FileProcessFactory:  # pylint: disable=too-few-public-methods
    SUPPORTED_EXTENSIONS = {
        'h5': ['.h5ad'],
        'csv': ['.csv']
    }

    @classmethod
    def get(cls, filename: str, file: BytesIO, **kwargs):
        if any([filename.endswith(ext) for ext in cls.SUPPORTED_EXTENSIONS['h5']]):
            return H5ADFileProcessor(file, **kwargs)
        if any([filename.endswith(ext) for ext in cls.SUPPORTED_EXTENSIONS['csv']]):
            return CSVFileProcessor(file, **kwargs)
        raise WrongExtension
