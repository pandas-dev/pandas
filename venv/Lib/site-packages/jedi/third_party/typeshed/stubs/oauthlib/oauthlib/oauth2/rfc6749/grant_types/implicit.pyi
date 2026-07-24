from _typeshed import Incomplete
from logging import Logger

from oauthlib.common import Request

from ..tokens import TokenBase
from .base import GrantTypeBase

log: Logger

class ImplicitGrant(GrantTypeBase):
    response_types: list[str]
    grant_allows_refresh_token: bool
    def create_authorization_response(
        self, request: Request, token_handler: TokenBase
    ) -> tuple[dict[str, str], str | None, int]: ...
    def create_token_response(self, request: Request, token_handler: TokenBase) -> tuple[dict[str, str], str | None, int]: ...
    def validate_authorization_request(self, request: Request) -> tuple[Incomplete, dict[str, Incomplete]]: ...
    def validate_token_request(self, request: Request) -> tuple[Incomplete, dict[str, Incomplete]]: ...
