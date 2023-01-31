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
        for extension in cls.EXTENSIONS_MAP.keys():
            if filename.endswith(extension):
                processor = cls.EXTENSIONS_MAP[extension]
                return processor(file, **kwargs)
        raise WrongExtension
