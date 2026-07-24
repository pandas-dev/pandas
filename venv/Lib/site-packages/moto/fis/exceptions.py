"""Exceptions raised by the fis service."""

from moto.core.exceptions import JsonRESTError


class FISException(JsonRESTError):
    pass


class ResourceNotFoundException(FISException):
    code = 404

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)
