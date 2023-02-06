from typing import Optional

from anndata import AnnData

from file_process.exceptions import NoColumnsError, ModelFileValidationVariablesError
from file_process.h5ad.schemas import SbioModelData


class H5ADValidator:
    def __init__(self, adata: AnnData, model_data: Optional[SbioModelData] = None):
        self.adata = adata
        self.model_data = model_data

    def __call__(self):
        self.validate()
        self.model_file_validation()

    def validate(self):
        target_names = list(self.adata.obs.columns)
        if not target_names:
            raise NoColumnsError

    def model_file_validation(self):
        if not self.model_data:
            return
        dataset_vars = set(self.adata.var.index)
        difference = self.model_data.var_names - dataset_vars
        if difference:
            raise ModelFileValidationVariablesError(difference)
