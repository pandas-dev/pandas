from _typeshed import Incomplete
from logging import Logger

from oauthlib.common import _HTTPMethod

from .base import BaseEndpoint

log: Logger

class AuthorizationEndpoint(BaseEndpoint):
    def __init__(self, default_response_type: str, default_token_type: str, response_types: dict[str, Incomplete]) -> None: ...
    @property
    def response_types(self) -> dict[str, Incomplete]: ...
    @property
    def default_response_type(self) -> str: ...
    @property
    def default_response_type_handler(self): ...
    @property
    def default_token_type(self): ...
    def create_authorization_response(
        self,
        uri: str,
        http_method: _HTTPMethod = "GET",
        body: str | None = None,
        headers: dict[str, str] | None = None,
        scopes=None,
        credentials: dict[str, Incomplete] | None = None,
    ): ...
    def validate_authorization_request(
        self, uri: str, http_method: _HTTPMethod = "GET", body: str | None = None, headers: dict[str, str] | None = None
    ): ...
