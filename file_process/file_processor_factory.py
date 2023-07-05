from io import BytesIO

from file_process.csv.csv_processor import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.h5ad.h5ad_processor import H5ADFileProcessor


class FileProcessFactory:  # pylint: disable=too-few-public-methods
    EXTENSIONS_MAP = {
        '.h5ad': H5ADFileProcessor,
        '.h5': H5ADFileProcessor,
        '.h5Seurat': H5ADFileProcessor,
        '.csv': CSVFileProcessor,
        '.txt': CSVFileProcessor,
        '.tsv': CSVFileProcessor,
    }

    @classmethod
    def get(cls, filename: str, file: BytesIO, **kwargs):
        for extension, processor_class in cls.EXTENSIONS_MAP.items():
            if filename.endswith(extension):
                kwargs['ext']=extension
                return processor_class(file, **kwargs)
        raise WrongExtension

    @classmethod
    def validate_extension(cls, filename: str):
        if not filename:
            raise WrongExtension
        for extension in cls.EXTENSIONS_MAP:
            if filename.endswith(extension):
                return
        raise WrongExtension
