from abc import ABC, abstractmethod
from typing import List
from io import BytesIO

import pandas as pd
from anndata import AnnData

from data_validation_preview.exceptions import ModelFileValidationError


class FileProcessorBase(ABC):
    @abstractmethod
    def read_file(self, file, **kwargs):
        raise NotImplemented

    @abstractmethod
    def process(self, file, model_csv_file: BytesIO = None, **kwargs) -> (List[str], List[dict], List[dict]):
        raise NotImplemented


class TabularFileProcessorBase(FileProcessorBase, ABC):
    def get_preview_data(self, adata):
        target_names = list(adata.obs.columns)
        var_preview = adata.var.head(n=10)
        obs_preview = adata.obs.head(n=10)
        return target_names, var_preview, obs_preview

    def model_file_validation(self, adata: AnnData, model_csv_file: BytesIO):
        reader = pd.read_csv(model_csv_file, sep=',', index_col=0)
        var_names = reader.index
        dataset_vars = list(adata.var.index)
        result = all(elem in dataset_vars for elem in var_names)
        if not result:
            raise ModelFileValidationError

    def process(self, file, model_csv_file: BytesIO = None, **kwargs) -> (List[str], List[dict], List[dict]):
        adata = self.read_file(file, **kwargs)
        if model_csv_file:
            self.model_file_validation(adata, model_csv_file)
        target_names, var_preview, obs_preview = self.get_preview_data(adata)
        return target_names, var_preview, obs_preview
