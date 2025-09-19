"""Exceptions raised by the connectcampaigns service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    """When a resource is not found."""

    code = 404

    def __init__(self, message: str) -> None:
        super().__init__("ResourceNotFoundException", message)


class ValidationException(JsonRESTError):
    """When validation fails on input parameters."""

    code = 400

    def __init__(self, message: str) -> None:
        super().__init__("ValidationException", message)
