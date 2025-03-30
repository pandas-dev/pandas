"""Exceptions raised by the emrcontainers service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, resource: str):
        super().__init__("ResourceNotFoundException", resource)
