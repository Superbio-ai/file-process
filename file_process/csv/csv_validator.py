import json
from io import BytesIO
from typing import Optional

from file_process.exceptions import NotAllTargetsError, NotSomeTargetsError, ModelFileValidationVariablesError, \
    CustomValidationException
from pandas import DataFrame

from file_process.csv.schemas import TabularValidationRules


class CSVValidator:
    def __init__(self, data_df: DataFrame, validation_rules: Optional[dict],
                 model_metadata_file: Optional[BytesIO] = None):
        self.data_df = data_df
        self.rules = TabularValidationRules(validation_rules) if validation_rules else None
        self.model_metadata_file = model_metadata_file

    def __call__(self):
        self.validate()
        self.model_file_validation()

    def validate(self):
        if not self.rules:
            return
        self._validate_column_names()
        self._validate_columns()

    def _validate_column_names(self):
        if self.rules.column_names_required:
            column_names = self.data_df.columns.values.tolist()
            for col in self.rules.columns:
                if col.required and col.name not in column_names:
                    all_required_columns = [col.name for col in self.rules.columns if col.required]
                    raise CustomValidationException(f'Missing {col.name} column in the file. '
                                                    f'List of required columns: [{all_required_columns}]')
            if not self.rules.accept_other_columns:
                allowed_column_names = [col.name for col in self.rules.columns]
                for col in column_names:
                    if col not in allowed_column_names:
                        raise CustomValidationException(f'Invalid column: {col}. '
                                                        f'The list of allowed column names: {allowed_column_names}')
        else:
            raise NotImplemented

    def _validate_columns(self):
        rules = {c.name: c for c in self.rules.columns}
        for name, data in self.data_df.iteritems():
            self._validate_column(name, data, rules.get(name))

    def _validate_column(self, name, data, rule):
        if rule.allowed_types:
            type_ = rule.allowed_types[0]
            try:
                data.astype(type_)
            except Exception as e:
                text = str(e)
                raise CustomValidationException(f'All values under {name} column must be one of the following types: '
                                                f'{rule.allowed_types}. {text.capitalize()}.')
        if not rule.allow_missings:
            if data.isna().sum():
                raise CustomValidationException(f'Column {name} has missings and it is not allowed.')
        if not rule.allow_duplicates:
            if data.duplicated().any():
                raise CustomValidationException(f'Column {name} has duplicates and it is not allowed.')
        if rule.min is not None:
            if data.le(rule.min).any():
                raise CustomValidationException(f'Min value in column {name} can be {rule.min}.')
        if rule.max is not None:
            if data.le(rule.max).any():
                raise CustomValidationException(f'Max value in column {name} can be {rule.max}.')
        if rule.allowed_values:
            if not data.eq(rule.allowed).all():
                raise CustomValidationException(f'For {name} column the list of allowed values is {rule.allowed_values}.')

    def model_file_validation(self):
        if not self.model_metadata_file:
            return
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
