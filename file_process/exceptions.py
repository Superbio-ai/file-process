from abc import ABC
from typing import Set


class BaseError(ABC, Exception):
    pass


class NotAllTargetsError(BaseError):
    def __init__(self, targets_missing: Set[str]):
        self.message = f'Validation check failed: new data does not contain all targets required by model. ' \
                       f'Missing targets: {targets_missing}'


class NotSomeTargetsError(BaseError):
    def __init__(self, model_targets: Set[str]):
        self.message = f'Validation check failed: new data does not contain any targets required by model. ' \
                       f'List of targets in a model: {model_targets}'


class ModelFileValidationVariablesError(BaseError):
    def __init__(self, variables_missing: Set[str]):
        self.message = f'Validation check failed: new data does not contain all variables (columns) required by ' \
                       f'model. Missing variables: {variables_missing}'


class WrongExtension(BaseError):
    message = 'Cannot process file: wrong extension.'


class DelimiterError(BaseError):
    message = 'Parsing error: try changing delimiter.'


class NoColumnsError(BaseError):
    message = 'No columns in file.'


class CustomValidationException(BaseError):
    def __init__(self, message: str):
        self.message = message


class NoXExpression(BaseError):
    message = 'The h5ad artifact does not contain expression data ".X".'


class DataIsNormalized(BaseError):
    message = 'Data should not be normalized.'


class DataIsNotFinite(BaseError):
    message = 'Data must be finite.'
