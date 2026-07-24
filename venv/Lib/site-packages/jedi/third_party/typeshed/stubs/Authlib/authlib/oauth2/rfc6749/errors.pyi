from _typeshed import Incomplete

from authlib.oauth2 import OAuth2Error as OAuth2Error

__all__ = [
    "OAuth2Error",
    "InsecureTransportError",
    "InvalidRequestError",
    "InvalidClientError",
    "UnauthorizedClientError",
    "InvalidGrantError",
    "UnsupportedResponseTypeError",
    "UnsupportedGrantTypeError",
    "InvalidScopeError",
    "AccessDeniedError",
    "MissingAuthorizationError",
    "UnsupportedTokenTypeError",
    "MissingCodeException",
    "MissingTokenException",
    "MissingTokenTypeException",
    "MismatchingStateException",
]

class InsecureTransportError(OAuth2Error):
    error: str
    description: str
    @classmethod
    def check(cls, uri) -> None: ...

class InvalidRequestError(OAuth2Error):
    error: str

class InvalidClientError(OAuth2Error):
    error: str
    status_code: int
    def get_headers(self): ...

class InvalidGrantError(OAuth2Error):
    error: str

class UnauthorizedClientError(OAuth2Error):
    error: str

class UnsupportedResponseTypeError(OAuth2Error):
    error: str
    response_type: Incomplete
    def __init__(
        self,
        response_type,
        description=None,
        uri=None,
        status_code=None,
        state=None,
        redirect_uri=None,
        redirect_fragment: bool = False,
        error=None,
    ) -> None: ...
    def get_error_description(self): ...

class UnsupportedGrantTypeError(OAuth2Error):
    error: str
    grant_type: Incomplete
    def __init__(self, grant_type) -> None: ...
    def get_error_description(self): ...

class InvalidScopeError(OAuth2Error):
    error: str
    description: str

class AccessDeniedError(OAuth2Error):
    error: str
    description: str

class ForbiddenError(OAuth2Error):
    status_code: int
    auth_type: Incomplete
    realm: Incomplete
    def __init__(self, auth_type=None, realm=None) -> None: ...
    def get_headers(self): ...

class MissingAuthorizationError(ForbiddenError):
    error: str
    description: str

class UnsupportedTokenTypeError(ForbiddenError):
    error: str

class MissingCodeException(OAuth2Error):
    error: str
    description: str

class MissingTokenException(OAuth2Error):
    error: str
    description: str

class MissingTokenTypeException(OAuth2Error):
    error: str
    description: str

class MismatchingStateException(OAuth2Error):
    error: str
    description: str
