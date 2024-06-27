import json
from typing import Any, List, Tuple

from moto.core.exceptions import JsonRESTError


class ManagedBlockchainClientError(JsonRESTError):
    def __init__(self, error_type: str, message: str):
        super().__init__(error_type=error_type, message=message)
        self.error_type = error_type
        self.message = message
        self.description = json.dumps({"message": self.message})

    def get_headers(self, *args: Any, **kwargs: Any) -> List[Tuple[str, str]]:  # pylint: disable=unused-argument
        return [
            ("Content-Type", "application/json"),
            ("x-amzn-ErrorType", self.error_type),
        ]

    def get_body(self, *args: Any, **kwargs: Any) -> str:  # pylint: disable=unused-argument
        return self.description


class BadRequestException(ManagedBlockchainClientError):
    def __init__(self, pretty_called_method: str, operation_error: str):
        super().__init__(
            "BadRequestException",
            f"An error occurred (BadRequestException) when calling the {pretty_called_method} operation: {operation_error}",
        )


class InvalidRequestException(ManagedBlockchainClientError):
    def __init__(self, pretty_called_method: str, operation_error: str):
        super().__init__(
            "InvalidRequestException",
            f"An error occurred (InvalidRequestException) when calling the {pretty_called_method} operation: {operation_error}",
        )


class ResourceNotFoundException(ManagedBlockchainClientError):
    def __init__(self, pretty_called_method: str, operation_error: str):
        self.code = 404
        super().__init__(
            "ResourceNotFoundException",
            f"An error occurred (ResourceNotFoundException) when calling the {pretty_called_method} operation: {operation_error}",
        )


class ResourceAlreadyExistsException(ManagedBlockchainClientError):
    def __init__(self, pretty_called_method: str, operation_error: str):
        self.code = 409
        super().__init__(
            "ResourceAlreadyExistsException",
            f"An error occurred (ResourceAlreadyExistsException) when calling the {pretty_called_method} operation: {operation_error}",
        )


class ResourceLimitExceededException(ManagedBlockchainClientError):
    def __init__(self, pretty_called_method: str, operation_error: str):
        self.code = 429
        super().__init__(
            "ResourceLimitExceededException",
            f"An error occurred (ResourceLimitExceededException) when calling the {pretty_called_method} operation: {operation_error}",
        )
