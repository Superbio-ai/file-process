from io import BytesIO
from typing import List

import anndata
import pandas as pd

from file_process.base import FileProcessorBase
from file_process.exceptions import ModelFileValidationError, NoColumnsError


class H5ADFileProcessor(FileProcessorBase):

    def __init__(self, file, **_):
        self.adata = anndata.read_h5ad(file)

    def get_targets(self):
        return list(self.adata.obs.columns)

    def get_observations(self):
        return self.adata.obs.head(n=10)

    def get_variables(self):
        return self.adata.var.head(n=10)

    def get_preview(self):
        target_names = self.get_targets()
        obs_preview = self.get_observations()
        var_preview = self.get_variables()
        return target_names, self.create_tabular_response(obs_preview), self.create_tabular_response(var_preview)

    def model_file_validation(self, model_metadata_file: BytesIO):
        reader = pd.read_csv(model_metadata_file, sep=',', index_col=0)
        var_names = reader.index
        dataset_vars = list(self.adata.var.index)
        result = all(elem in dataset_vars for elem in var_names)
        if not result:
            raise ModelFileValidationError

    def validate(self):
        target_names = list(self.adata.obs.columns)
        if not target_names:
            raise NoColumnsError

    @staticmethod
    def create_tabular_response(data_df: pd.DataFrame) -> List[dict]:
        if data_df is None:
            return []
        data_df = data_df.astype(object)
        rows = data_df.round(2).where(pd.notnull(data_df), None).to_dict(orient='records')
        indices = list(data_df.index)
        if len(rows) != len(indices):
            return rows
        for index, value in enumerate(indices):
            rows[index]['Feature Name'] = value
        return rows
