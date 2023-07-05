import pytest

from file_process.exceptions import ModelFileValidationVariablesError, NoColumnsError, NoXExpression, \
    DataIsNormalized, DataIsNotFinite
from file_process.h5ad.h5ad_processor import H5ADFileProcessor
from file_process.h5ad.h5ad_validator import H5ADValidator
from file_process.h5ad.schemas import SbioModelDataForH5ad
from tests.test_file_process import H5AD_INPUT_FILES_PATH, get_remote_file_obj


class TestH5ADValidator:
    path = f'{H5AD_INPUT_FILES_PATH}/heart_sample.h5ad'
    valid_model_path = f'{H5AD_INPUT_FILES_PATH}/heart_sample_model.csv'
    invalid_model_path = f'{H5AD_INPUT_FILES_PATH}/heart_sample_invalid_model.csv'

    def test_h5ad_model_file_validation(self):
        file_bytes_io = get_remote_file_obj(self.path)
        model_data = SbioModelDataForH5ad(get_remote_file_obj(self.valid_model_path))
        ext = self.path.split('.')[-1]
        processor = H5ADFileProcessor(file_bytes_io, ext)
        validator = H5ADValidator(processor.adata, model_data)
        _ = validator.model_file_validation()

    def test_h5ad_model_file_validation_invalid_model(self):
        file_bytes_io = get_remote_file_obj(self.path)
        model_data = SbioModelDataForH5ad(get_remote_file_obj(self.invalid_model_path))
        ext = self.path.split('.')[-1]
        processor = H5ADFileProcessor(file_bytes_io, ext)
        validator = H5ADValidator(processor.adata, model_data)
        with pytest.raises(ModelFileValidationVariablesError):
            _ = validator.model_file_validation()

    def test_validate_target_names_present_invalid(self):
        file_bytes_io = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/heart_atlas_no_columns.h5ad')
        processor = H5ADFileProcessor(file_bytes_io, ext='h5ad')
        validator = H5ADValidator(processor.adata)
        with pytest.raises(NoColumnsError):
            validator._validate_target_names_present()

    def test_validate_target_names_present_valid(self):
        file_bytes_io = get_remote_file_obj(self.path)
        ext = self.path.split('.')[-1]
        processor = H5ADFileProcessor(file_bytes_io, ext)
        validator = H5ADValidator(processor.adata)
        validator._validate_target_names_present()

    def test_validate_check_x_invalid(self):
        file_bytes_io = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/heart_atlas_no_x.h5ad')
        processor = H5ADFileProcessor(file_bytes_io, ext='h5ad')
        validator = H5ADValidator(processor.adata)
        with pytest.raises(NoXExpression):
            validator._check_x()

    def test_validate_check_x_valid(self):
        file_bytes_io = get_remote_file_obj(self.path)
        ext = self.path.split('.')[-1]
        processor = H5ADFileProcessor(file_bytes_io, ext)
        validator = H5ADValidator(processor.adata)
        validator._check_x()

    def test_validate_check_not_normed(self):
        file_bytes_io = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/heart_sample_normalized.h5ad')
        processor = H5ADFileProcessor(file_bytes_io, ext='h5ad')
        validator = H5ADValidator(processor.adata)
        with pytest.raises(DataIsNormalized):
            validator._check_not_normed()

    def test_validate_check_normed_valid(self):
        file_bytes_io = get_remote_file_obj(self.path)
        ext = self.path.split('.')[-1]
        processor = H5ADFileProcessor(file_bytes_io, ext)
        validator = H5ADValidator(processor.adata)
        validator._check_not_normed()

    def test_validate_check_finite_invalid(self):
        file_bytes_io = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/heart_atlas_nan_inf.h5ad')
        processor = H5ADFileProcessor(file_bytes_io, ext='h5ad')
        validator = H5ADValidator(processor.adata)
        with pytest.raises(DataIsNotFinite):
            validator._check_finite()

    def test_validate_check_finite_valid(self):
        file_bytes_io = get_remote_file_obj(self.path)
        ext = self.path.split('.')[-1]
        processor = H5ADFileProcessor(file_bytes_io, ext)
        validator = H5ADValidator(processor.adata)
        validator._check_finite()
