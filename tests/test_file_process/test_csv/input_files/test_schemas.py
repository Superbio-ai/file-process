import pytest
from pydantic import ValidationError

from file_process.csv.schemas import TabularValidationRules


class TestSchemas:
    def test_validations(self):
        valid_config = {
            'columns': {
                'columnsList': [
                    {'name': 'Samantha', 'required': True, 'allowedTypes': ['str'],
                     'allowedValues': ['mac', 'and', 'cheese'], 'allowMissings': True, 'allowDuplicates': True},
                    {'name': 'Karrie', 'required': True, 'allowedTypes': ['int'],
                     'allowedValues': [1, 2, 3], 'allowMissings': True, 'allowDuplicates': True}
                ],
                'columnNamesRequired': True,
                'allowOtherColumns': True,
                'preserveOrder': False
            }
        }

        validation_rules = TabularValidationRules(valid_config)
        validation_rules.validate_self()

    invalid_columns = [
        {'name': 'test', 'allowedTypes': ['str'], 'allowedValues': [1, 2, 3]},  # wrong type of allowed values
        # {'allowedTypes': ['str'], 'allowedValues': ['strings']},  # no name
        # {'name': '', 'allowedTypes': ['str'], 'allowedValues': ['strings']},  # empty name
        # {'name': 'test', 'allowedTypes': ['str', 'int'], 'allowedValues': ['strings']},  # 2 allowed types
        {'name': 'test', 'allowedTypes': ['int'], 'min': 0.1},  # wrong type of min
        {'name': 'test', 'allowedTypes': ['int'], 'max': 0.1},  # wrong type of max
        {'name': 'test', 'allowedTypes': ['int'], 'min': 10, 'max': 0},  # min bugger than max
        {'name': 'test', 'allowedTypes': ['int'], 'max': 0, 'min': 10},  # min bigger than max - different fields order
    ]

    @pytest.mark.parametrize('invalid_column', invalid_columns)
    def test_validations_invalid_columns(self, invalid_column):
        invalid_config = {
            'columns': {
                'columnsList': [invalid_column],
                'preserveOrder': False
            }
        }
        validation_rules = TabularValidationRules(invalid_config)
        errors = validation_rules.validate_self()
        assert errors
