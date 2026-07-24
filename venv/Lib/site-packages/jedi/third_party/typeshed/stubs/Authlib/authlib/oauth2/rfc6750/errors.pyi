from _typeshed import Incomplete

from authlib.oauth2 import OAuth2Error

__all__ = ["InvalidTokenError", "InsufficientScopeError"]

class InvalidTokenError(OAuth2Error):
    error: str
    description: str
    status_code: int
    realm: Incomplete
    extra_attributes: Incomplete
    def __init__(self, description=None, uri=None, status_code=None, state=None, realm=None, **extra_attributes) -> None: ...
    def get_headers(self) -> list[tuple[str, str]]: ...

class InsufficientScopeError(OAuth2Error):
    error: str
    description: str
    status_code: int
