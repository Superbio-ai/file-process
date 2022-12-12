from io import BytesIO
from typing import List

import anndata
import pandas as pd
from numpy import nan

from file_process.exceptions import ModelFileValidationError, NoColumnsError

from file_process.base import FileProcessorBase


class H5ADFileProcessor(FileProcessorBase):

    def __init__(self, file, **kwargs):
        self.data = anndata.read_h5ad(file)

    def get_targets(self) -> List[str]:
        return list(self.data.obs.columns)

    def get_obs(self):
        return self.data.obs.head(n=10)

    def get_var(self):
        return self.data.var.head(n=10)

    def get_preview(self) -> (List[str], List[dict], List[dict]):
        target_names = self.get_targets()
        try:
            obs = self.get_obs()
            obs_preview = self._create_tabular_response(obs)
        except Exception as exc:
            obs_preview = None

        try:
            var = self.get_var()
            var_preview = self._create_tabular_response(var)
        except Exception as exc:
            var_preview = None

        return target_names, obs_preview, var_preview

    def model_file_validation(self, model_metadata_file: BytesIO):
        reader = pd.read_csv(model_metadata_file, sep=',', index_col=0)
        var_names = reader.index
        dataset_vars = list(self.data.var.index)
        result = all(elem in dataset_vars for elem in var_names)
        if not result:
            raise ModelFileValidationError

    def validate(self):
        target_names = self.get_targets()
        if not target_names:
            raise NoColumnsError

    @staticmethod
    def _create_tabular_response(data_df: pd.DataFrame) -> List[dict]:
        if data_df is None:
            return []
        rows = data_df.round(2).replace({nan: None}).to_dict(orient='records')
        indices = list(data_df.index)
        for index, value in enumerate(indices):
            rows[index]['Feature Name'] = value
        return rows
