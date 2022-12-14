from io import BytesIO


H5AD_INPUT_FILES_PATH = 'tests/test_file_process/input_files/h5ad'
CSV_INPUT_FILES_PATH = 'tests/test_file_process/input_files/tabular_csv'


def get_remote_file_obj(path: str):
    with open(path, 'rb') as file:
        file_obj = BytesIO(file.read())
        return file_obj
