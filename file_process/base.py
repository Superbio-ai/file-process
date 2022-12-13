from abc import ABC, abstractmethod
from io import BytesIO
from typing import List, Optional


class FileProcessorBase(ABC):

    @abstractmethod
    def get_preview(self) -> (List[str], List[dict], Optional[List[dict]]):
        pass

    @abstractmethod
    def model_file_validation(self, model_file: BytesIO) -> None:
        pass

    @abstractmethod
    def validate(self):
        pass
