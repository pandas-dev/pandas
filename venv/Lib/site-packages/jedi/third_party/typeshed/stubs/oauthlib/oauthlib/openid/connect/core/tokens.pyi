from collections.abc import Callable

from oauthlib.common import Request
from oauthlib.oauth2.rfc6749.tokens import TokenBase as TokenBase

from .request_validator import RequestValidator

class JWTToken(TokenBase):
    __slots__ = ("request_validator", "token_generator", "refresh_token_generator", "expires_in")
    request_validator: RequestValidator
    token_generator: Callable[[Request], str] | Callable[[Request, bool], str]
    refresh_token_generator: Callable[[Request], str] | Callable[[Request, bool], str]
    expires_in: int | Callable[[Request], int]
    def __init__(
        self,
        request_validator: RequestValidator | None = None,
        token_generator: Callable[[Request], str] | None = None,
        expires_in: int | Callable[[Request], int] | None = None,
        refresh_token_generator: Callable[[Request], str] | None = None,
    ) -> None: ...
    def create_token(self, request: Request, refresh_token: bool = False): ...
    def validate_request(self, request: Request): ...
    def estimate_type(self, request: Request): ...
