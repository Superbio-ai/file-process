import json
from io import BytesIO
from typing import Optional

from file_process.exceptions import NotAllTargetsError, NotSomeTargetsError, ModelFileValidationVariablesError
from pandas import DataFrame

from file_process.csv.schemas import CSVValidationRules


class CSVValidator:
    def __int__(self, data_df: DataFrame, validation_rules: CSVValidationRules,
                model_metadata_file: Optional[BytesIO] = None):
        self.data_df = data_df
        self.validation_rules = validation_rules
        self.model_metadata_file = model_metadata_file

    def __call__(self):
        self.validate()
        self.model_file_validation()

    def validate(self):
        if self.validation_rules.column_names_required:
            # TODO check that columns are named
            # TODO make a validation based on column names
            pass
        else:
            # TODO make a validation based on indexes
            pass

    def model_file_validation(self):
        reader = json.load(self.model_metadata_file)
        var_names = set(reader['columns'])
        target_names = set(reader['targets'])
        metadata = reader.get('metadata', {})
        dataset_vars = set(self.data_df.columns)

        all_targets = metadata.get('require_all_targets', 'all')
        if all_targets == 'all':
            difference = target_names - dataset_vars
            if difference:
                raise NotAllTargetsError(difference)
        elif all_targets == 'some':
            are_targets_valid = not target_names or any(elem in dataset_vars for elem in target_names)
            if not are_targets_valid:
                raise NotSomeTargetsError(target_names)
        dataset_diff = dataset_vars - target_names
        var_names_diff = var_names - target_names
        difference = var_names_diff - dataset_diff
        if difference:
            raise ModelFileValidationVariablesError(difference)

