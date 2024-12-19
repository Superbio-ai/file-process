from typing import Optional

from pandas import DataFrame

from file_process.exceptions import NotAllTargetsError, NotSomeTargetsError, ModelFileValidationVariablesError, \
    CustomValidationException
from file_process.csv.schemas import TabularValidationRules, SbioModelDataForCsv, ColumnValidationRule


class CSVValidator:
    def __init__(self, data_df: DataFrame, validation_rules: Optional[dict],
                 model_data: Optional[SbioModelDataForCsv] = None):
        self.data_df = data_df
        self.rules = TabularValidationRules(validation_rules)
        self.model_data = model_data

    def __call__(self):
        self.validate()
        self.model_file_validation()

    def validate(self):
        if not self.rules:
            return
        if self.rules.column_names_required:
            self._validate_column_names()
            self._validate_columns()

    def _validate_column_names(self):
        column_names = self.data_df.columns.values.tolist()
        if self.rules.preserve_order:
            # all columns are considered as required with this filter
            required_columns = [col.name for col in self.rules.columns]
            if required_columns != column_names:
                raise CustomValidationException(f'This is the list of allowed columns in the allowed order: '
                                                f'[{required_columns}]')
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

    def _validate_columns(self):
        rules = {c.name: c for c in self.rules.columns}
        for name, data in self.data_df.items():
            if name in rules:
                self._validate_column(str(name), data, rules[name])

    def _validate_column(self, name: str, data, rule: ColumnValidationRule):
        if rule.allowed_types:
            # TODO add support for a list of types
            type_ = rule.allowed_types[0]
            try:
                data.astype(type_)
            except Exception as e:
                text = str(e)
                raise CustomValidationException(f'All values under {name} column must be one of the following '
                                                f'types: {rule.allowed_types}. {text.capitalize()}.') from e
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
            if data.ge(rule.max).any():
                raise CustomValidationException(f'Max value in column {name} can be {rule.max}.')
        if rule.allowed_values:
            if not data.isin(rule.allowed_values).all():
                raise CustomValidationException(f'For {name} column the list of allowed values is '
                                                f'{rule.allowed_values}.')

    def model_file_validation(self):
        if not self.model_data:
            return

        model_target_names = self.model_data.target_names
        dataset_vars = set(self.data_df.columns)
        all_targets = self.model_data.metadata.get('require_all_targets', 'all')

        if all_targets == 'all':
            difference = model_target_names - dataset_vars
            if difference:
                raise NotAllTargetsError(difference)
        elif all_targets == 'some':
            are_targets_valid = not model_target_names or any(elem in dataset_vars for elem in model_target_names)
            if not are_targets_valid:
                raise NotSomeTargetsError(model_target_names)
        dataset_diff = dataset_vars - model_target_names
        var_names_diff = self.model_data.var_names - model_target_names
        difference = var_names_diff - dataset_diff
        if difference:
            raise ModelFileValidationVariablesError(difference)
