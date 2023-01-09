from abc import ABC, abstractmethod
from io import BytesIO
from typing import List


class FileProcessorBase(ABC):
    @abstractmethod
    def validate(self):
        pass

    @abstractmethod
    def model_file_validation(self, model_metadata_file: BytesIO):
        pass

    @abstractmethod
    def get_preview(self) -> (List[str], List[dict], List[dict]):
        pass
