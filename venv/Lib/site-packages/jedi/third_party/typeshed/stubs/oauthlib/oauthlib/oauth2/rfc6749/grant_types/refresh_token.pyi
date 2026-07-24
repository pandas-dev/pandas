from collections.abc import Iterable
from logging import Logger

from oauthlib.common import Request

from ..request_validator import RequestValidator
from ..tokens import TokenBase
from .base import GrantTypeBase, _AuthValidator, _TokenValidator

log: Logger

class RefreshTokenGrant(GrantTypeBase):
    def __init__(
        self,
        request_validator: RequestValidator | None = None,
        issue_new_refresh_tokens: bool = True,
        *,
        post_auth: Iterable[_AuthValidator] | None = None,
        post_token: Iterable[_TokenValidator] | None = None,
        pre_auth: Iterable[_AuthValidator] | None = None,
        pre_token: Iterable[_TokenValidator] | None = None,
        **kwargs,
    ) -> None: ...
    def create_token_response(self, request: Request, token_handler: TokenBase) -> tuple[dict[str, str], str, int | None]: ...
    def validate_token_request(self, request: Request) -> None: ...
