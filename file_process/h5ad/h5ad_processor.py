from io import BytesIO
from typing import List, Optional

import anndata
import pandas as pd
import numpy as np
import h5py

from file_process.base import FileProcessorBase
from file_process.constants import PREVIEW_ROWS_COUNT
from file_process.h5ad.h5ad_validator import H5ADValidator
from file_process.h5ad.schemas import SbioModelDataForH5ad


class H5GenericProcessor(FileProcessorBase):

    def validate(self, model_metadata_file: Optional[BytesIO] = None, _: Optional[dict] = None):
        model_data = None
        if model_metadata_file:
            model_data = SbioModelDataForH5ad(model_metadata_file)
        validator = H5ADValidator(self.adata, model_data, enable_warnings=False)
        validator()

    def get_targets(self):
        return list(self.adata.obs.columns)

    def get_observations(self, rows_number: Optional[int] = None):
        return self.adata.obs.head(n=rows_number)

    def get_variables(self, rows_number: Optional[int] = None):
        return self.adata.var.head(n=rows_number)

    def get_preview(self):
        target_names = self.get_targets()
        obs_preview = self.get_observations(PREVIEW_ROWS_COUNT)
        var_preview = self.get_variables(PREVIEW_ROWS_COUNT)
        return target_names, self.create_tabular_response(obs_preview), self.create_tabular_response(var_preview)

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


class H5FileProcessor(H5GenericProcessor):

    def __init__(self, file, **_):
        data_mat = h5py.File(file, 'r')
        main_df = np.array(data_mat['X'])
        if 'Y' in data_mat:
            obs = np.array(data_mat['Y'])
        else:
            obs = None
        data_mat.close()

        # preprocessing scRNA-seq read counts matrix
        adata = anndata.AnnData(main_df)
        if obs is not None:
            adata.obs['Group'] = obs

        self.adata = adata


class H5ADFileProcessor(H5GenericProcessor):

    def __init__(self, file, **_):
        self.adata = anndata.read_h5ad(file)
