from typing import Any

from moto.core.exceptions import JsonRESTError


class MacieException(JsonRESTError):
    pass


class ValidationException(MacieException):
    code = 400

    def __init__(self, message: str):
        super().__init__("ValidationException", message)

    def get_headers(self, *args: Any, **kwargs: Any) -> list[tuple[str, str]]:
        return [("x-amzn-ErrorType", "ValidationException")]


class AccessDeniedException(MacieException):
    code = 403

    def __init__(self, message: str):
        super().__init__("AccessDeniedException", message)

    def get_headers(self, *args: Any, **kwargs: Any) -> list[tuple[str, str]]:
        return [("x-amzn-ErrorType", "AccessDeniedException")]


class ResourceNotFoundException(MacieException):
    code = 404

    def __init__(self, message: str):
        super().__init__("ResourceNotFoundException", message)

    def get_headers(self, *args: Any, **kwargs: Any) -> list[tuple[str, str]]:
        return [("x-amzn-ErrorType", "ResourceNotFoundException")]
