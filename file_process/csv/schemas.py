class ColumnValidationRule:
    type_map = {
        'str': str,
        'int': int,
        'float': float
    }

    def __init__(self, validation_rules: dict):
        self.name = validation_rules.get('name')
        # self.allowed_types = [self.type_map[type_] for type_ in validation_rules.get('allowedTypes')]
        self.allowed_types = validation_rules.get('typesAllowed')
        self.required = validation_rules.get('required')
        self.allow_missings = validation_rules.get('allowMissingss')
        self.allow_duplicates = validation_rules.get('allowDuplicates')
        self.min = validation_rules.get('min')
        self.max = validation_rules.get('max')
        self.allowed_values = validation_rules.get('allowedValues')


class TabularValidationRules:
    def __init__(self, validation_rules: dict):
        column_rules = validation_rules.get('columns')
        self.preserve_order = column_rules.get('preserveOrder')
        self.column_names_required = column_rules.get('columnNamesRequired')
        self.accept_other_columns = column_rules.get('allowOtherColumns')
        self.columns = [ColumnValidationRule(column_data) for column_data in column_rules.get('columnsList', [])]
