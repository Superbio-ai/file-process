from abc import ABC, abstractmethod
from typing import List

import pandas as pd


class FileProcessorBase(ABC):

    @abstractmethod
    def get_preview(self):
        pass

    @abstractmethod
    def model_file_validation(self, model_file):
        pass

