from abc import ABC, abstractmethod


class BaseError(ABC, Exception):
    @property
    @abstractmethod
    def message(self):
        pass


class ModelFileValidationError(BaseError):
    message = 'Validation check failed: new data does not contain all fields required by model.'


class ModelFileValidationTargetsError(BaseError):
    message = 'Validation check failed: new data does not contain all targets required by model.'


class ModelFileValidationVariablesError(BaseError):
    message = 'Validation check failed: new data does not contain all variables (columns) required by model.'


class WrongExtension(BaseError):
    message = 'Cannot process file: wrong extension.'


class DelimiterError(BaseError):
    message = 'Parsing error: try changing delimiter.'


class NoColumnsError(BaseError):
    message = 'No columns in file.'
