"""Exceptions raised by the emrserverless service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, resource: str, resource_type: str = "Application"):
        super().__init__(
            "ResourceNotFoundException", f"{resource_type} {resource} does not exist"
        )


class ValidationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)


class AccessDeniedException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__("AccessDeniedException", message)
