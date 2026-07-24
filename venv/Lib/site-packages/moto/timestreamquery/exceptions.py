"""Exceptions raised by the timestreamquery service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFound(JsonRESTError):
    def __init__(self, arn: str):
        super().__init__(
            error_type="ResourceNotFoundException",
            message=f"The resource with arn {arn} does not exist.",
        )
