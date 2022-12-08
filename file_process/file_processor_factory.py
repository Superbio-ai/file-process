from file_process.csv import CSVFileProcessor
from file_process.exceptions import WrongExtension
from file_process.h5 import H5FileProcessor


class FileProcessFactory:  # pylint: disable=too-few-public-methods
    @classmethod
    def get(cls, filename: str, file):
        if filename.endswith('.h5') or filename.endswith('.h5ad'):
            h5_processor = H5FileProcessor()
            try:
                # File extension alone is not enough to determine the read method:
                # sometimes h5-10x files may have .h5 extension.
                _ = h5_processor.read_file(file)
                return h5_processor
            except TypeError as exc:
                raise WrongExtension from exc
        if filename.endswith('.csv'):
            return CSVFileProcessor()
        raise WrongExtension
