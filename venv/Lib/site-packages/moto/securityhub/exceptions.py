"""Exceptions raised by the securityhub service."""

from moto.core.exceptions import ServiceException


class SecurityHubClientError(ServiceException):
    pass


class _InvalidOperationException(SecurityHubClientError):
    def __init__(self, error_type: str, op: str, msg: str):
        super().__init__(
            error_type,
            f"An error occurred ({error_type}) when calling the {op} operation: {msg}",
        )


class InvalidInputException(_InvalidOperationException):
    def __init__(self, op: str, msg: str):
        super().__init__("InvalidInputException", op, msg)
