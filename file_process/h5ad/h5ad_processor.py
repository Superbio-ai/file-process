import os
import tempfile
from io import BytesIO
from typing import List, Optional

import anndata
import pandas as pd

from file_process.base import FileProcessorBase
from file_process.constants import PREVIEW_ROWS_COUNT
from file_process.h5ad.h5ad_validator import H5ADValidator
from file_process.h5ad.schemas import SbioModelDataForH5ad


class H5ADFileProcessor(FileProcessorBase):

    def __init__(self, file, **_):
        if isinstance(file, BytesIO):
            with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as temp_file:
                temp_file.write(file.getvalue())
                temp_file_path = temp_file.name
        elif isinstance(file, str):
            temp_file_path = file
        else:
            raise TypeError(f"Unsupported file type: {type(file)}")

        self.adata = anndata.read_h5ad(temp_file_path)

        if isinstance(file, BytesIO):
            os.remove(temp_file_path)

    def validate(self, model_metadata_file: Optional[BytesIO] = None, _: Optional[dict] = None):
        model_data = None
        if model_metadata_file:
            model_data = SbioModelDataForH5ad(model_metadata_file)
        validator = H5ADValidator(self.adata, model_data, enable_warnings=False)
        validator()

    def get_obs_targets(self):
        return list(self.adata.obs.columns)

    def get_var_targets(self):
        return list(self.adata.var.columns)

    def get_observations(self, rows_number: Optional[int] = None):
        return self.adata.obs.head(n=rows_number)

    def get_variables(self, rows_number: Optional[int] = None):
        return self.adata.var.head(n=rows_number)

    def get_preview(self):
        var_target_names = self.get_var_targets()
        obs_target_names = self.get_obs_targets()
        obs_preview = self.get_observations(PREVIEW_ROWS_COUNT)
        var_preview = self.get_variables(PREVIEW_ROWS_COUNT)
        return var_target_names, obs_target_names, self.create_tabular_response(
            obs_preview), self.create_tabular_response(var_preview), None

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
