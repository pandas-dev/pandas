"""Exceptions raised by the connect service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    code = 404

    def __init__(self, message: str) -> None:
        super().__init__("ResourceNotFoundException", message)


class InvalidParameterException(JsonRESTError):
    code = 400

    def __init__(self, message: str) -> None:
        super().__init__("InvalidParameterException", message)
