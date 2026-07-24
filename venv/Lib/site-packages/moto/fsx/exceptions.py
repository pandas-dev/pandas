"""Exceptions raised by the fsx service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    def __init__(self, msg: str):
        super().__init__("ResourceNotFoundException", f"{msg}")
