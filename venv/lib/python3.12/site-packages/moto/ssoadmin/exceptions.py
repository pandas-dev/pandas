"""Exceptions raised by the ssoadmin service."""

from moto.core.exceptions import JsonRESTError


class ResourceNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, message: str = "") -> None:
        super().__init__(
            error_type="ResourceNotFoundException",
            message=message,
        )


class ConflictException(JsonRESTError):
    code = 400

    def __init__(self, message: str = "") -> None:
        super().__init__(
            error_type="ConflictException",
            message=message,
        )


class ServiceQuotaExceededException(JsonRESTError):
    code = 400

    def __init__(self, message: str = "") -> None:
        super().__init__(
            error_type="ServiceQuotaExceededException",
            message=message,
        )
