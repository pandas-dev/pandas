"""Exceptions raised by the workspaces service."""

from moto.core.exceptions import JsonRESTError


class ValidationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class InvalidParameterValuesException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidParameterValuesException", message)


class ResourceAlreadyExistsException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceAlreadyExistsException", message)


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)
