import json
from io import BytesIO
from typing import List

import pandas as pd

NAME_STR = 'name'
ALLOWED_TYPES_STR = 'allowedTypes'
ALLOWED_VALUES_STR = 'allowedValues'
REQUIRED_STR = 'required'
MIN_STR = 'min'
MAX_STR = 'max'
PRESERVE_ORDER_STR = 'preserveOrder'
COLUMN_NAMES_REQUIRED = 'columnNamesRequired'
ALLOW_OTHER_COLUMNS_STR = 'allowOtherColumns'
COLUMNS_LIST_STR = 'columnsList'


class ValidationRuleError:
    default_path = 'columns'

    def __init__(self, keys: List[str], message: str):
        self.keys = [f'{self.default_path}.{key}' for key in keys]
        self.message = message

    def to_dict(self):
        return self.__dict__


class ColumnValidationRule:

    def __init__(self, validation_rules: dict):
        self._validation_rules = validation_rules
        self.name = validation_rules.get(NAME_STR)
        self.allowed_types = validation_rules.get(ALLOWED_TYPES_STR)
        self.required = validation_rules.get(REQUIRED_STR, True)
        self.allow_missings = validation_rules.get('allowMissings', True)
        self.allow_duplicates = validation_rules.get('allowDuplicates', True)
        self.min = validation_rules.get(MIN_STR)
        self.max = validation_rules.get(MAX_STR)
        self.allowed_values = validation_rules.get(ALLOWED_VALUES_STR)

    def validate_self(self, index: int):  # pylint: disable=too-many-branches
        errors = []
        for key in self._validation_rules.keys():
            if key not in [NAME_STR, ALLOWED_TYPES_STR, REQUIRED_STR, 'allowMissings', 'allowDuplicates', MIN_STR,
                           MAX_STR, ALLOWED_VALUES_STR]:
                errors.append(ValidationRuleError([f'{COLUMNS_LIST_STR}[{index}].{key}'], 'Invalid key'))
        if not self.name:
            errors.append(ValidationRuleError(
                [f'{COLUMNS_LIST_STR}[{index}].{NAME_STR}'],
                'Missing field'
            ))
        if self.allowed_types and len(self.allowed_types) > 1:
            errors.append(ValidationRuleError(
                [f'{COLUMNS_LIST_STR}[{index}].{ALLOWED_TYPES_STR}'],
                'We only support one allowed type at the moment (notice that float includes int).'
            ))
        if self.min is not None and self.max is not None:
            try:
                if self.min > self.max:
                    errors.append(ValidationRuleError(
                        [f'{COLUMNS_LIST_STR}[{index}].{MIN_STR}', f'{COLUMNS_LIST_STR}[{index}].{MAX_STR}'],
                        'Min cannot be bigger than max.'
                    ))
            except TypeError as e:
                errors.append(ValidationRuleError(
                    [f'{COLUMNS_LIST_STR}[{index}].{MIN_STR}', f'{COLUMNS_LIST_STR}[{index}].{MAX_STR}'],
                    f'Impossible to do comparison between min and max. Error: {e}'
                ))
        if self.min is not None:
            try:
                self._validate_type(self.min)
            except ValueError:
                errors.append(ValidationRuleError([f'{COLUMNS_LIST_STR}[{index}].{MIN_STR}'],
                                                  'Min value must be one of allowed types.'))
        if self.max is not None:
            try:
                self._validate_type(self.max)
            except ValueError:
                errors.append(ValidationRuleError([f'{COLUMNS_LIST_STR}[{index}].{MAX_STR}'],
                                                  'Max value must be one of allowed types.'))
        if self.allowed_types and len(self.allowed_types) == 1 and self.allowed_values:
            try:
                self._validate_type(self.allowed_values)
            except ValueError:
                errors.append(ValidationRuleError(
                    [f'{COLUMNS_LIST_STR}[{index}].{ALLOWED_VALUES_STR}',
                     f'{COLUMNS_LIST_STR}[{index}].{ALLOWED_TYPES_STR}'],
                    'All allowed values must be one of allowed type.'
                ))
        return errors

    def _validate_type(self, value):
        if not self.allowed_types:
            return
        fixed_value = [value] if not isinstance(value, list) else value
        values_df = pd.DataFrame(fixed_value)
        values_df.astype(self.allowed_types[0], errors='raise')


class TabularValidationRules:
    def __init__(self, validation_rules: dict):
        if not validation_rules:
            validation_rules = {}
        column_rules = validation_rules.get('columns', {})
        self._column_rules = column_rules
        self.preserve_order = column_rules.get(PRESERVE_ORDER_STR)
        self.column_names_required = column_rules.get(COLUMN_NAMES_REQUIRED, True)
        self.accept_other_columns = column_rules.get(ALLOW_OTHER_COLUMNS_STR, True)
        self.columns = [ColumnValidationRule(column_data) for column_data in column_rules.get('columnsList', [])]

    def validate_self(self) -> List[ValidationRuleError]:
        errors = []
        for key in self._column_rules.keys():
            if key not in [PRESERVE_ORDER_STR, COLUMN_NAMES_REQUIRED, ALLOW_OTHER_COLUMNS_STR, 'columnsList']:
                errors.append(ValidationRuleError([key], 'Invalid key'))
        if not self.column_names_required:
            errors.append(ValidationRuleError(
                [COLUMN_NAMES_REQUIRED],
                'Field must be always true. Validation by index is not implemented yet. '
            ))
        if self.accept_other_columns and self.preserve_order:
            errors.append(ValidationRuleError(
                [PRESERVE_ORDER_STR, ALLOW_OTHER_COLUMNS_STR],
                f'If you want to set preserveOrder to true, you need to explicitly set {ALLOW_OTHER_COLUMNS_STR} '
                f'to false. When we validate {PRESERVE_ORDER_STR} we get the list of all allowed columns '
                f'and compare it with the list of all columns in the file, so if there will be other columns, '
                f'the validation will fail.'
            ))

        for index, column in enumerate(self.columns):
            if self.preserve_order and not column.required:
                errors.append(ValidationRuleError(
                    [PRESERVE_ORDER_STR, f'{COLUMNS_LIST_STR}[{index}].{REQUIRED_STR}'],
                    f'If {PRESERVE_ORDER_STR} is true, then all columns must be required.'
                ))
            errors += column.validate_self(index)
        return errors


class SbioModelDataForCsv:
    # TODO make one class for csv and h5ad
    def __init__(self, model_metadata_file: BytesIO):
        reader = json.load(model_metadata_file)
        self.var_names = set(reader['columns'])
        self.target_names = set(reader['targets'])
        self.metadata = reader.get('metadata', {})
