import json
from io import BytesIO


class ColumnValidationRule:

    def __init__(self, validation_rules: dict):
        self.name = validation_rules.get('name')
        self.allowed_types = validation_rules.get('typesAllowed')
        self.required = validation_rules.get('required', True)
        self.allow_missings = validation_rules.get('allowMissingss', True)
        self.allow_duplicates = validation_rules.get('allowDuplicates', True)
        self.min = validation_rules.get('min')
        self.max = validation_rules.get('max')
        self.allowed_values = validation_rules.get('allowedValues')


class TabularValidationRules:
    def __init__(self, validation_rules: dict):
        if not validation_rules:
            validation_rules = {}
        column_rules = validation_rules.get('columns', {})
        self.preserve_order = column_rules.get('preserveOrder')
        self.column_names_required = column_rules.get('columnNamesRequired', True)
        self.accept_other_columns = column_rules.get('allowOtherColumns', True)
        self.columns = [ColumnValidationRule(column_data) for column_data in column_rules.get('columnsList', [])]


class SbioModelDataForCsv:
    # TODO make one class for csv and h5ad
    def __init__(self, model_metadata_file: BytesIO):
        reader = json.load(model_metadata_file)
        self.var_names = set(reader['columns'])
        self.target_names = set(reader['targets'])
        self.metadata = reader.get('metadata', {})
