from collections.abc import Iterable
from logging import Logger

from oauthlib.common import Request
from oauthlib.oauth2.rfc6749.grant_types.base import _AuthValidator, _TokenValidator
from oauthlib.oauth2.rfc6749.grant_types.implicit import ImplicitGrant as OAuth2ImplicitGrant
from oauthlib.oauth2.rfc6749.request_validator import RequestValidator as OAuth2RequestValidator

from .base import GrantTypeBase

log: Logger

class ImplicitGrant(GrantTypeBase):
    proxy_target: OAuth2ImplicitGrant
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
    def openid_authorization_validator(self, request: Request): ...
