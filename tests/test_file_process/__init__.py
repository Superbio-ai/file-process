from io import BytesIO


H5AD_INPUT_FILES_PATH = 'tests/test_file_process/test_h5ad/input_files'
CSV_INPUT_FILES_PATH = 'tests/test_file_process/test_csv/input_files'


def get_remote_file_obj(path: str):
    with open(path, 'rb') as file:
        file_obj = BytesIO(file.read())
        return file_obj
