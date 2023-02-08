class ColumnValidationRule:
    def __init__(self, validation_rules: dict):
        self.name = validation_rules.get('name')
        self.allowed_types = validation_rules.get('allowedTypes')
        self.required = validation_rules.get('required')
        self.allow_missings = validation_rules.get('allowMissingss')
        self.allow_duplicates = validation_rules.get('allowDuplicates')
        self.min = validation_rules.get('min')
        self.max = validation_rules.get('max')
        self.allowed_values = validation_rules.get('allowedValues')


class CSVValidationRules:
    def __init__(self, validation_rules: dict):
        self.preserve_order = validation_rules.get('preserveOrder')
        self.column_names_required = validation_rules.get('columnNamesRequired')
        self.accept_other_columns = validation_rules.get('allowOtherColumns')
        self.columns = [ColumnValidationRule(column_data) for column_data in validation_rules.get('columns', [])]
