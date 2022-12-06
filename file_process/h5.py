import anndata

from file_process.base import TabularFileProcessorBase


class H5FileProcessor(TabularFileProcessorBase):

    def read_file(self, file, **kwargs):
        return anndata.read_h5ad(file)