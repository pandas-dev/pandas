from _typeshed import Incomplete
from logging import Logger

from oauthlib.common import Request
from oauthlib.oauth2.rfc6749.request_validator import RequestValidator as OAuth2RequestValidator

log: Logger

class Dispatcher:
    default_grant: Incomplete | None
    oidc_grant: Incomplete | None

class AuthorizationCodeGrantDispatcher(Dispatcher):
    default_grant: Incomplete | None
    oidc_grant: Incomplete | None
    def __init__(self, default_grant=None, oidc_grant=None) -> None: ...
    def create_authorization_response(self, request: Request, token_handler): ...
    def validate_authorization_request(self, request: Request): ...

class ImplicitTokenGrantDispatcher(Dispatcher):
    default_grant: Incomplete | None
    oidc_grant: Incomplete | None
    def __init__(self, default_grant=None, oidc_grant=None) -> None: ...
    def create_authorization_response(self, request: Request, token_handler): ...
    def validate_authorization_request(self, request: Request): ...

class AuthorizationTokenGrantDispatcher(Dispatcher):
    default_grant: Incomplete | None
    oidc_grant: Incomplete | None
    request_validator: OAuth2RequestValidator
    def __init__(self, request_validator: OAuth2RequestValidator, default_grant=None, oidc_grant=None) -> None: ...
    def create_token_response(self, request: Request, token_handler): ...
