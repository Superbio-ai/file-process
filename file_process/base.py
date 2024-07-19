from abc import ABC, abstractmethod
from io import BytesIO
from typing import Optional, List


class FileProcessorBase(ABC):
    @abstractmethod
    def validate(self, model_metadata_file: Optional[BytesIO] = None, validation_rules: Optional[dict] = None):
        pass

    @abstractmethod
    def get_preview(self) -> (List[str], List[str], List[dict], List[dict], str):
        pass
