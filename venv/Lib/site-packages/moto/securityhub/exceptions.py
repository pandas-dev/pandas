"""Exceptions raised by the securityhub service."""

from moto.core.exceptions import JsonRESTError


class SecurityHubClientError(JsonRESTError):
    code = 400


class _InvalidOperationException(SecurityHubClientError):
    def __init__(self, error_type: str, op: str, msg: str):
        super().__init__(
            error_type,
            "An error occurred (%s) when calling the %s operation: %s"
            % (error_type, op, msg),
        )


class InvalidInputException(_InvalidOperationException):
    def __init__(self, op: str, msg: str):
        super().__init__("InvalidInputException", op, msg)
