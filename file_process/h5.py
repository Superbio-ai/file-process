from io import BytesIO

import anndata
import scanpy as sc
from werkzeug.datastructures import FileStorage

from file_process.base import TabularFileProcessorBase


class H5FileProcessor(TabularFileProcessorBase):

    def read_file(self, file, **kwargs):
        return anndata.read_h5ad(file)


class H510xFileProcessor(TabularFileProcessorBase):

    def read_file(self, file, **kwargs):
        if isinstance(file, FileStorage):  # TODO try to get rid of it
            file = file.read()
            # file = BytesIO(file)
        return sc.read_10x_h5(file)
