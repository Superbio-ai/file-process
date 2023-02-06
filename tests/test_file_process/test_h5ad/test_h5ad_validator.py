import pytest

from file_process.exceptions import ModelFileValidationVariablesError, NoColumnsError
from file_process.h5ad.h5ad_processor import H5ADFileProcessor
from file_process.h5ad.h5ad_validator import H5ADValidator
from file_process.h5ad.schemas import SbioModelData
from tests.test_file_process import H5AD_INPUT_FILES_PATH, get_remote_file_obj


class TestH5ADValidator:
    path = f'{H5AD_INPUT_FILES_PATH}/heart_sample.h5ad'
    valid_model_path = f'{H5AD_INPUT_FILES_PATH}/heart_sample_model.csv'
    invalid_model_path = f'{H5AD_INPUT_FILES_PATH}/heart_sample_invalid.csv'

    def test_h5ad_model_file_validation(self):
        file_bytes_io = get_remote_file_obj(self.path)
        model_data = SbioModelData(get_remote_file_obj(self.valid_model_path))
        processor = H5ADFileProcessor(file_bytes_io)
        validator = H5ADValidator(processor.adata, model_data)
        _ = validator.model_file_validation()

    def test_h5ad_model_file_validation_invalid_model(self):
        file_bytes_io = get_remote_file_obj(self.path)
        model_data = SbioModelData(get_remote_file_obj(self.invalid_model_path))
        processor = H5ADFileProcessor(file_bytes_io)
        validator = H5ADValidator(processor.adata, model_data)
        with pytest.raises(ModelFileValidationVariablesError):
            _ = validator.model_file_validation()

    def test_validate_no_columns(self):
        file_bytes_io = get_remote_file_obj(f'{H5AD_INPUT_FILES_PATH}/pbmc3k_raw.h5ad')
        processor = H5ADFileProcessor(file_bytes_io)
        validator = H5ADValidator(processor.adata)
        with pytest.raises(NoColumnsError):
            validator.validate()
