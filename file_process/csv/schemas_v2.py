import json
from io import BytesIO
from typing import List, Optional

import pandas as pd
from pydantic import BaseModel, Field, validator, root_validator  # TODO add pydantic to requirements


class ColumnValidationRule(BaseModel):
    name: str = Field(min_length=1)
    required: Optional[bool] = Field(default=True)
    allowed_types: Optional[List[str]] = Field(alias='allowedTypes', max_items=1)
    allowed_values: Optional[List[str]] = Field(alias='allowedValues')
    allow_missings: Optional[bool] = Field(alias='allowMissings', default=True)
    allow_duplicates: Optional[bool] = Field(alias='allowDuplicates', default=True)
    min: Optional[float]
    max: Optional[float]

    class Config:
        allow_population_by_field_name = True

    @validator('allowed_values')
    def name_length(cls, allowed_values, values):
        try:
            cls._validate_type(allowed_values, values)
        except Exception as e:
            raise ValueError('All allowed values must be one of the allowed types.')
        return allowed_values

    def _validate_type(cls, value, values):
        if not value:
            return
        allowed_type = values['allowed_types'][0] if values.get('allowed_types') else None
        values = [value] if not isinstance(value, list) else value
        values_df = pd.DataFrame(values)
        values_df.astype(allowed_type)

    @validator('max')
    def valid_max(cls, value, values):
        try:
            cls._validate_type(value, values)
        except Exception as e:
            raise ValueError(f'Max has type that is different from allowed types.')
        return value

    @validator('min')
    def valid_min(cls, value, values):
        try:
            cls._validate_type(value, values)
        except Exception as e:
            raise ValueError(f'Mib has type that is different from allowed types.')
        return value

    @root_validator()
    def validate_all_fields_at_the_same_time(cls, values):
        if values.get('min') is None or values.get('max') is None:
            return values
        assert values.get('min') <= values.get('max'), 'Min cannot be bigger than max.'
        return values


class TabularValidationRules(BaseModel):
    columns: Optional[List[ColumnValidationRule]] = Field(alias='columnsList')
    column_names_required: Optional[bool] = Field(default=True, alias='columnNamesRequired')
    accept_other_columns: Optional[bool] = Field(default=True, alias='allowOtherColumns')
    preserve_order: Optional[bool] = Field(default=False, alias='preserveOrder')

    @validator('column_names_required')
    def names_required_validation(cls, value):
        if not value:
            raise ValueError('Field must be always true. Validation by index is not implemented yet.')
        return value

    @validator('preserve_order')
    def preserve_order_validation(cls, preserve_order, values):
        if not preserve_order:
            return preserve_order
        if values['accept_other_columns']:
            raise ValueError('Unable to set both preserveOrder and allowOtherColumns to true '
                             'because it does not make sense')
        for column in values['columns']:
            if not column.required:
                raise ValueError('If preserveOrder is true, then all columns must be required.')


class SbioModelDataForCsv:
    # TODO make one class for csv and h5ad
    def __init__(self, model_metadata_file: BytesIO):
        reader = json.load(model_metadata_file)
        self.var_names = set(reader['columns'])
        self.target_names = set(reader['targets'])
        self.metadata = reader.get('metadata', {})
