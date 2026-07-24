from ..base import OAuth2Error

__all__ = ["InvalidRequestUriError", "InvalidRequestObjectError", "RequestNotSupportedError", "RequestUriNotSupportedError"]

class InvalidRequestUriError(OAuth2Error):
    error: str
    description: str
    status_code: int

class InvalidRequestObjectError(OAuth2Error):
    error: str
    description: str
    status_code: int

class RequestNotSupportedError(OAuth2Error):
    error: str
    description: str
    status_code: int

class RequestUriNotSupportedError(OAuth2Error):
    error: str
    description: str
    status_code: int
