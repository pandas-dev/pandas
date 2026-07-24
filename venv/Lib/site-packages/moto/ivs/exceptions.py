"""Exceptions raised by the ivs service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    code = 404

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)
