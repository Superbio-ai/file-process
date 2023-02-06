from typing import Optional

import numpy as np
from anndata import AnnData
from scanpy.get import _get_obs_rep
from scipy.sparse import issparse

from file_process.exceptions import NoColumnsError, ModelFileValidationVariablesError
from file_process.h5ad.schemas import SbioModelData
from file_process.logger import logger


class H5ADValidator:
    def __init__(self, adata: AnnData, model_data: Optional[SbioModelData] = None, enable_warnings: bool = True):
        self.adata = adata
        self.model_data = model_data
        self.enable_warnings = enable_warnings

    def __call__(self):
        self.validate()
        self.model_file_validation()

    def validate(self):
        self._validate_target_names_present()
        self._check_x()
        self._check_normed()
        self._check_finite()

        if self.enable_warnings:
            warnings = []
            warnings += self._check_structure_warnings()
            warnings += self._validate_encoding_version()
            if warnings:
                logger.info(f"Warnings: {warnings}")

    def model_file_validation(self):
        if not self.model_data:
            return
        dataset_vars = set(self.adata.var.index)
        difference = self.model_data.var_names - dataset_vars
        if difference:
            raise ModelFileValidationVariablesError(difference)

    def _validate_target_names_present(self):
        target_names = list(self.adata.obs.columns)
        if not target_names:
            raise NoColumnsError

    def _check_x(self):
        if not hasattr(self.adata, 'X'):
            raise Exception('The h5ad artifact does not contain expression data ".X".')

    def _check_structure_warnings(self):
        warnings = []
        if not hasattr(self.adata, 'obs'):
            warnings.append(
                'The h5ad artifact does not contain observation information ".obs".'
            )

        if not hasattr(self.adata, 'var'):
            warnings.append(
                'The h5ad artifact does not contain variable information ".var".'
            )

        if not hasattr(self.adata, 'obsm'):
            warnings.append(
                'The h5ad artifact does not contain experiment design information ".obsm".'
            )

        if not hasattr(self.adata, 'uns'):
            warnings.append(
                'The h5ad artifact does not contain schema information ".uns".'
            )
        return warnings

    def _check_normed(self, obs_key: Optional[str] = None):
        data = _get_obs_rep(self.adata, layer=obs_key)
        diff_sum = np.array(data != data[0]).sum()
        if diff_sum == 0:
            raise Exception('Data cannot be normalized.')

    def _check_finite(self, obs_key: Optional[str] = None):
        data = _get_obs_rep(self.adata, layer=obs_key)
        if issparse(data):
            is_finite = np.isfinite(data.A).all()
        else:
            is_finite = np.isfinite(data).all()
        if is_finite:
            raise Exception('The data is finite.')

    def _validate_encoding_version(self):
        pass
        # encoding_dict = dict(f.attrs)
        # encoding_version = encoding_dict.get("encoding-version")
        # if encoding_version != "0.1.0":
        #     "The h5ad artifact was generated with an AnnData version different from 0.8.0."
