"""Exceptions raised by the shield service."""

from moto.core.exceptions import JsonRESTError


class ResourceAlreadyExistsException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceAlreadyExistsException", message)


class InvalidResourceException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidResourceException", message)


class InvalidParameterException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("InvalidParameterException", message)


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)


class ValidationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)
