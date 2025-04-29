"""Exceptions raised by the quicksight service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("ResourceNotFoundException", msg)


class InvalidParameterValueException(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("InvalidParameterValueException", msg)


class ValidationException(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("ValidationException", msg)


class SchemaException(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("SchemaException", msg)


class ParamValidationError(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("ParamValidationError", msg)
