from io import BytesIO

from file_process.csv.csv_processor import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.h5ad.h5ad_processor import H5ADFileProcessor
from file_process.txt.txt_processor import TxtFileProcessor


class FileProcessFactory:  # pylint: disable=too-few-public-methods
    EXTENSIONS_MAP = {
        '.h5ad': H5ADFileProcessor,
        '.csv': CSVFileProcessor,
        '.tsv': CSVFileProcessor,
        '.txt': TxtFileProcessor,
        '.pdb': TxtFileProcessor,
        '.fasta': TxtFileProcessor,
    }

    @classmethod
    def get(cls, filename: str, file: BytesIO, **kwargs):
        for extension, processor_class in cls.EXTENSIONS_MAP.items():
            if filename.endswith(extension):
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
