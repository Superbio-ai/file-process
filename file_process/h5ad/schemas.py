from io import BytesIO

import pandas as pd


class SbioModelDataForH5ad:
    def __init__(self, model_metadata_file: BytesIO):
        reader = pd.read_csv(model_metadata_file, sep=',', index_col=0)
        self.var_names = set(reader.index)
