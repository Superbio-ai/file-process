import pytest
from pytest_mock import MockerFixture
from werkzeug.datastructures import FileStorage
import pandas as pd
from numpy import nan
import json
from numpy import number as numpy_number

from data_validation_preview.csv import CSVFileProcessor
from data_validation_preview.exceptions import ModelFileValidationError
from tests.test_data_validation_preview import INPUT_FILES_PATH


class TestCSVFileProcessor:
    main_path = f'{INPUT_FILES_PATH}/tabular_csv/'
    all_configs_path = f'{main_path}all_configs.json'
    validity_path = f'{main_path}validity.json'
    mock_configs_path = f'{main_path}mock_configs/'
    mock_data_path = f'{main_path}mock_data/'
    original_data_path = f'{main_path}original_data.csv'
    preview_data_path = f'{main_path}preview_data/original_data_preview.json'
    preview_variables_path = f'{main_path}preview_data/original_variables_preview.json'

    def read_json(self, _p):
        with open(_p) as json_file:
            data = json.load(json_file)
        return data
    
    def write_json(self, dictionary, name):
        json_object = json.dumps(dictionary, indent=4)
        with open(name, "w") as outfile:
            outfile.write(json_object)

    def create_valid_configs(self):
        '''Used to create test data and configs'''
        data = pd.read_csv(self.original_data_path)
        path = self.mock_configs_path
        name1 = f'{self.mock_data_path}valid_new_data.csv'
        name2 = f'{self.mock_data_path}sometimes_valid_new_data.csv'
        name3 = f'{self.mock_data_path}invalid_new_data.csv'   
        data.sample(n=data.shape[1],axis='columns',replace=False).to_csv(name1, index=False)
        data[list(data.columns[0:3]) + [data.columns[4]]].to_csv(name2, index=False)
        data[data.columns[1:]].to_csv(name3, index=False)
        all_configs = []
        validity = {
            name1:[],
            name2:[],
            name3:[]
        }
        if path != None:
            path = path if path[-1] == '/' else path + '/'
        else:
            path = ''
        # ##Preview Data
        preview_data = data.head(min(10, data.shape[0]))
        original_columns = list(preview_data.columns)
        numeric_columns = preview_data.select_dtypes(include=numpy_number).columns
        preview_data = preview_data.replace({nan: None})
        preview_data[numeric_columns] = preview_data[numeric_columns].round(2)  
        preview_data = preview_data.to_dict(orient='records')
        self.write_json(preview_data, self.preview_data_path)
        self.write_json(original_columns, self.preview_variables_path)
        # ## Config Files
        # Test 1
        nname = path + 'valid_train_supervised_data_config_1.json'
        valid_train_supervised_data_config_1 = {
            'columns':  list(data.columns),
            'targets': list([data.columns[3]])
        }
        self.write_json(valid_train_supervised_data_config_1,nname)
        all_configs.append(nname)
        validity[name1].append(nname)
        # Test 2
        nname = path + 'valid_train_supervised_data_config_2.json'
        valid_train_supervised_data_config_2 = {
            'columns': list(data.columns),
            'targets': list([data.columns[3]]),
            'metadata':{
                'require_all_targets':True
            }
        }
        self.write_json(valid_train_supervised_data_config_2,nname)
        all_configs.append(nname)
        validity[name1].append(nname)
        # Test 3
        nname = path +  'valid_train_supervised_data_config_3.json'
        valid_train_supervised_data_config_3 = {
            'columns':  list(data.columns),
            'targets': list(data.columns[2:4])
        }
        self.write_json(valid_train_supervised_data_config_3,nname)
        all_configs.append(nname)
        validity[name1].append(nname)
        # Test 4
        nname = path + 'valid_train_supervised_data_config_4.json'
        valid_train_supervised_data_config_4 = {
            'columns': list(data.columns),
            'targets': list(data.columns[2:4]),
            'metadata':{
                'require_all_targets':False
            }
        }
        self.write_json(valid_train_supervised_data_config_4,nname)
        all_configs.append(nname)
        validity[name1].append(nname)
        validity[name2].append(nname)

        # Test 5
        nname = path + 'valid_train_supervised_data_config_5.json'
        targs = list(data.columns[2:4])
        valid_train_supervised_data_config_5 = {
            'columns': [c for c in data.columns if c not in targs],
            'targets': targs,
            'metadata':{
                'require_all_targets':True
            }
        }
        self.write_json(valid_train_supervised_data_config_5,nname)
        all_configs.append(nname)
        validity[name1].append(nname)
        # Test 6
        nname = path + 'valid_usupervised_data_config_1.json'
        valid_usupervised_data_config_1 = {
            'columns': list(data.columns),
            'targets': [],
            'metadata':{
                'require_all_targets':False
            }
        }
        self.write_json(valid_usupervised_data_config_1,nname)
        all_configs.append(nname)
        validity[name1].append(nname)
        # Test 7
        nname = path + 'valid_usupervised_data_config_2.json'
        valid_usupervised_data_config_2 = {
            'columns': list(data.columns),
            'targets': [],
            'metadata':{
                'require_all_targets':True
            }
        }
        self.write_json(valid_usupervised_data_config_2,nname)
        all_configs.append(nname)
        validity[name1].append(nname)
        
        self.write_json(all_configs, self.all_configs_path)
        self.write_json(validity, self.validity_path)

    def _get_file_and_remote_file_obj(self, path: str):
        with open(path, 'rb') as file:
            file_obj = file.read()
            return file, file_obj

    def test_read_local_file(self):
        with open(self.original_data_path, 'rb') as file:
            file_obj = FileStorage(file)
            res = CSVFileProcessor().read_file(file_obj)
            assert isinstance(res, pd.DataFrame)

    def test_read_remote_file(self):
        file, file_obj = self._get_file_and_remote_file_obj(self.original_data_path)
        res = CSVFileProcessor().read_file(file_obj)
        assert isinstance(res, pd.DataFrame)
        file.close()

    def test_process(self):
        original_data_preview = self.read_json(self.preview_data_path)
        original_var_names = self.read_json(self.preview_variables_path)
        file, file_obj = self._get_file_and_remote_file_obj(self.original_data_path)
        var_names, obs_preview = CSVFileProcessor().process(file_obj)
        assert obs_preview == original_data_preview
        assert all(var in var_names for var in original_var_names)
        file.close()

    def tst_model_file_validation_with_csv(self, new_data_path, data_metadata_path):
        metadata_file, metadata_file_obj = self._get_file_and_remote_file_obj(data_metadata_path)
        test_file, test_file_obj = self._get_file_and_remote_file_obj(new_data_path)
        df = CSVFileProcessor().read_file(test_file_obj)
        is_valid = CSVFileProcessor().model_file_validation(df, metadata_file_obj)
        metadata_file.close()
        test_file.close()

    def tst_model_file_validation_with_csv_error(self, new_data_path, data_metadata_path):
        metadata_file, metadata_file_obj = self._get_file_and_remote_file_obj(data_metadata_path)
        test_file, test_file_obj = self._get_file_and_remote_file_obj(new_data_path)
        df = CSVFileProcessor().read_file(test_file_obj)
        with pytest.raises(ModelFileValidationError):
            is_valid = CSVFileProcessor().model_file_validation(df, metadata_file_obj)
        metadata_file.close()
        test_file.close()

    def test_csv_tabular_metadata_validation(self, mocker: MockerFixture):
        validity_dict = self.read_json(self.validity_path)
        all_configs = self.read_json(self.all_configs_path)  
        for dat in validity_dict:
            for met in all_configs:
                if met in validity_dict[dat]:
                    self.tst_model_file_validation_with_csv(dat, met, mocker)
                else:
                    self.tst_model_file_validation_with_csv_error(dat, met, mocker)
