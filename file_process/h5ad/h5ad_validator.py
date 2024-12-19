from typing import Optional, Union

import numpy as np
from anndata import AnnData
from scipy.sparse import issparse, spmatrix
import pandas as pd

from file_process.exceptions import NoColumnsError, ModelFileValidationVariablesError, NoXExpression, \
    DataIsNormalized, DataIsNotFinite
from file_process.h5ad.schemas import SbioModelDataForH5ad
from file_process.logger import logger


def _get_obs_rep(
        adata: AnnData,
        *,
        use_raw: bool = False,
        layer: Optional[str] = None,
        obsm: Optional[str] = None,
        obsp: Optional[str] = None,
) -> Union[np.ndarray, spmatrix, pd.DataFrame, None]:
    """
    Choose array aligned with obs annotation.
    """
    # https://github.com/scverse/scanpy/issues/1546
    if not isinstance(use_raw, bool):
        raise TypeError(f"use_raw expected to be bool, was {type(use_raw)}.")

    is_layer = layer is not None
    is_raw = use_raw is not False
    is_obsm = obsm is not None
    is_obsp = obsp is not None
    choices_made = sum((is_layer, is_raw, is_obsm, is_obsp))
    assert choices_made in {0, 1}
    if choices_made == 0:
        return adata.X
    if is_layer:
        return adata.layers[layer]
    if use_raw:
        return adata.raw.X
    if is_obsm:
        return adata.obsm[obsm]
    if is_obsp:
        return adata.obsp[obsp]
    raise AssertionError(
        "That was unexpected. Please report this bug at:\n\n\t"
        "https://github.com/scverse/scanpy/issues"
    )


class H5ADValidator:
    def __init__(self, adata: AnnData, model_data: Optional[SbioModelDataForH5ad] = None, enable_warnings: bool = True):
        self.adata = adata
        self.model_data = model_data
        self.enable_warnings = enable_warnings

    def __call__(self):
        self.validate()
        self.model_file_validation()

    def validate(self):
        self._validate_target_names_present()
        self._check_x()
        self._check_not_normed()
        self._check_finite()

        if self.enable_warnings:
            warnings = []
            warnings += self._get_structure_warnings()
            warnings += self._get_encoding_version_warnings()
            if warnings:
                logger.info("Warnings: %s", warnings)

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
        if not hasattr(self.adata, 'X') or self.adata.X is None:
            raise NoXExpression

    def _get_structure_warnings(self):
        warnings = []
        if not hasattr(self.adata, 'obs'):
            warnings.append('The h5ad artifact does not contain observation information ".obs".')

        if not hasattr(self.adata, 'var'):
            warnings.append('The h5ad artifact does not contain variable information ".var".')

        if not hasattr(self.adata, 'obsm'):
            warnings.append('The h5ad artifact does not contain experiment design information ".obsm".')

        if not hasattr(self.adata, 'uns'):
            warnings.append('The h5ad artifact does not contain schema information ".uns".')
        return warnings

    def _check_not_normed(self, obs_key: Optional[str] = None):
        data = _get_obs_rep(self.adata, layer=obs_key)
        if issparse(data):
            diff_sum = len(np.unique(np.round(data.toarray().sum(axis=1), 0), axis=0))
        else:
            diff_sum = len(np.unique(np.round(data.sum(axis=1), 0), axis=0))
        if diff_sum <= 2:
            raise DataIsNormalized

    def _check_finite(self, obs_key: Optional[str] = None):
        data = _get_obs_rep(self.adata, layer=obs_key)
        if issparse(data):
            is_finite = np.isfinite(data.toarray()).all()
        else:
            is_finite = np.isfinite(data).all()
        if not is_finite:
            raise DataIsNotFinite

    def _get_encoding_version_warnings(self):
        # TODO check
        return ["The h5ad artifact was generated with an AnnData version different from 0.8.0."]
