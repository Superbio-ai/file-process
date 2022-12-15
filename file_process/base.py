from abc import ABC, abstractmethod
from typing import List, Optional
from io import BytesIO

import pandas as pd
from anndata import AnnData
from pandas import notnull

from file_process.exceptions import ModelFileValidationError, NoColumnsError


class FileProcessorBase(ABC):
    @abstractmethod
    def read_file(self, file, **kwargs):
        raise NotImplemented

    @abstractmethod
    def process(self, file, model_metadata_file: BytesIO = None, **kwargs) \
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

    def model_file_validation(self, adata: AnnData, model_metadata_file: BytesIO):
        reader = pd.read_csv(model_metadata_file, sep=',', index_col=0)
        var_names = reader.index
        dataset_vars = list(adata.var.index)
        result = all(elem in dataset_vars for elem in var_names)
        if not result:
            raise ModelFileValidationError

    def validate(self, adata: AnnData):
        target_names = list(adata.obs.columns)
        if not target_names:
            raise NoColumnsError

    def process(self, file, model_metadata_file: BytesIO = None, **kwargs) -> (List[str], pd.DataFrame, pd.DataFrame):
        adata = self.read_file(file, **kwargs)
        self.validate(adata)
        if model_metadata_file:
            self.model_file_validation(adata, model_metadata_file)
        target_names, var_preview, obs_preview = self.get_preview_data(adata)
        return target_names, var_preview, obs_preview

    def create_tabular_response(self, data_df: pd.DataFrame) -> List[dict]:
        if data_df is None:
            return []
        data_df = data_df.astype(object)
        rows = data_df.round(2).where(notnull(data_df), None).to_dict(orient='records')
        indices = list(data_df.index)
        if len(rows) != len(indices):
            return rows
        for index, value in enumerate(indices):
            rows[index]['Feature Name'] = value
        return rows
