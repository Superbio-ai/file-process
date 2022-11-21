import pytest

from data_validation_preview.exceptions import WrongExtension
from data_validation_preview.file_processor_factory import FileProcessFactory
from data_validation_preview.h5 import H5FileProcessor


class TestFileProcessFactory:
    def test_get_h5_processor(self):
        res = FileProcessFactory.get('heart_atlas.h5ad')
        assert isinstance(res, H5FileProcessor)

    def test_get_h5ad_processor(self):
        res = FileProcessFactory.get('heart_atlas.h5')
        assert isinstance(res, H5FileProcessor)

    def test_get_not_implemented_processor(self):
        with pytest.raises(WrongExtension):
            _ = FileProcessFactory.get('heart_atlas.doc')
