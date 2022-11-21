import anndata

from data_validation_preview.base import TabularFileProcessorBase


class H5FileProcessor(TabularFileProcessorBase):

    def read_file(self, file, **kwargs):
        return anndata.read_h5ad(file)
