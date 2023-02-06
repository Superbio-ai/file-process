from abc import ABC, abstractmethod
from typing import List

from file_process.h5ad.schemas import SbioModelData


class FileProcessorBase(ABC):
    @abstractmethod
    def validate(self, model_data: SbioModelData):
        pass

    @abstractmethod
    def get_preview(self) -> (List[str], List[dict], List[dict]):
        pass
