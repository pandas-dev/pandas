from _typeshed import Incomplete
from logging import Logger

from oauthlib.common import Request, _HTTPMethod

from .base import BaseEndpoint

log: Logger

class ResourceEndpoint(BaseEndpoint):
    def __init__(self, default_token: str, token_types: dict[str, Incomplete]) -> None: ...
    @property
    def default_token(self) -> str: ...
    @property
    def default_token_type_handler(self): ...
    @property
    def tokens(self) -> dict[str, Incomplete]: ...
    def verify_request(
        self,
        uri: str,
        http_method: _HTTPMethod = "GET",
        body: str | None = None,
        headers: dict[str, str] | None = None,
        scopes=None,
    ) -> tuple[bool, Request]: ...
    def find_token_type(self, request: Request): ...
