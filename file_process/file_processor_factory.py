from io import BytesIO

from file_process.csv import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.h5 import H5ADFileProcessor


class FileProcessFactory:  # pylint: disable=too-few-public-methods
    EXTENSIONS_MAP = {
        '.h5ad': H5ADFileProcessor,
        '.csv': CSVFileProcessor
    }

    @classmethod
    def get(cls, filename: str, file: BytesIO, **kwargs):
        for extension, processor_class in cls.EXTENSIONS_MAP.items():
            if filename.endswith(extension):
                return processor_class(file, **kwargs)
        raise WrongExtension
