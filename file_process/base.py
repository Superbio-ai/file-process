from abc import ABC, abstractmethod
from typing import List, Optional
from io import BytesIO

import pandas as pd
from anndata import AnnData
from numpy import nan

from file_process.exceptions import ModelFileValidationError


class FileProcessorBase(ABC):
    @abstractmethod
    def read_file(self, file, **kwargs):
        raise NotImplemented

    @abstractmethod
    def process(self, file, model_csv_file: BytesIO = None, **kwargs) \
            -> (List[str], Optional[pd.DataFrame], Optional[pd.DataFrame]):
        raise NotImplemented

    @abstractmethod
    def create_tabular_response(self, data_df: pd.DataFrame) -> List[dict]:
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

    def process(self, file, model_csv_file: BytesIO = None, **kwargs) -> (List[str], pd.DataFrame, pd.DataFrame):
        adata = self.read_file(file, **kwargs)
        if model_csv_file:
            self.model_file_validation(adata, model_csv_file)
        target_names, var_preview, obs_preview = self.get_preview_data(adata)
        return target_names, var_preview, obs_preview

    def create_tabular_response(self, data_df: pd.DataFrame) -> List[dict]:
        rows = data_df.round(2).replace({nan: None}).to_dict(orient='records')
        indices = list(data_df.index)
        for index, value in enumerate(indices):
            rows[index]['Feature Name'] = value
        return rows
