from abc import ABC, abstractmethod
from typing import Set


class BaseError(ABC, Exception):
    pass


class ModelFileValidationError(BaseError):
    message = 'Validation check failed: new data does not contain all fields required by model.'


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
        self.message = f'Validation check failed: new data does not contain any variables (columns) required by model.' \
                       f'Missing variables: {variables_missing}'


class WrongExtension(BaseError):
    message = 'Cannot process file: wrong extension.'


class DelimiterError(BaseError):
    message = 'Parsing error: try changing delimiter.'


class NoColumnsError(BaseError):
    message = 'No columns in file.'
