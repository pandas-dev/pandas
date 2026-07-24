from collections.abc import Iterable
from logging import Logger

from oauthlib.common import Request
from oauthlib.oauth2.rfc6749.grant_types.base import _AuthValidator, _TokenValidator
from oauthlib.oauth2.rfc6749.grant_types.refresh_token import RefreshTokenGrant as OAuth2RefreshTokenGrant
from oauthlib.oauth2.rfc6749.request_validator import RequestValidator as OAuth2RequestValidator

from .base import GrantTypeBase

log: Logger

class RefreshTokenGrant(GrantTypeBase):
    proxy_target: OAuth2RefreshTokenGrant
    def __init__(
        self,
        request_validator: OAuth2RequestValidator | None = None,
        *,
        post_auth: Iterable[_AuthValidator] | None = None,
        post_token: Iterable[_TokenValidator] | None = None,
        pre_auth: Iterable[_AuthValidator] | None = None,
        pre_token: Iterable[_TokenValidator] | None = None,
        **kwargs,
    ) -> None: ...
    def add_id_token(self, token, token_handler, request: Request): ...  # type: ignore[override]
