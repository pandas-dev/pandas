from _typeshed import Incomplete
from collections.abc import Collection
from typing_extensions import TypeAlias

from authlib.oauth2 import OAuth2Request
from authlib.oauth2.rfc6749 import ClientMixin

from ..hooks import Hookable

_ServerResponse: TypeAlias = tuple[int, str, list[tuple[str, str]]]

class BaseGrant(Hookable):
    TOKEN_ENDPOINT_AUTH_METHODS: Collection[str]
    GRANT_TYPE: str | None
    TOKEN_RESPONSE_HEADER: Collection[tuple[str, str]]
    prompt: Incomplete
    redirect_uri: Incomplete
    request: OAuth2Request
    server: Incomplete
    def __init__(self, request: OAuth2Request, server) -> None: ...
    @property
    def client(self): ...
    def generate_token(
        self,
        user=None,
        scope: str | None = None,
        grant_type: str | None = None,
        expires_in: int | None = None,
        include_refresh_token: bool = True,
    ) -> dict[str, str | int]: ...
    def authenticate_token_endpoint_client(self) -> ClientMixin: ...
    def save_token(self, token): ...
    def validate_requested_scope(self) -> None: ...

class TokenEndpointMixin:
    TOKEN_ENDPOINT_HTTP_METHODS: Incomplete
    GRANT_TYPE: Incomplete
    @classmethod
    def check_token_endpoint(cls, request: OAuth2Request) -> bool: ...
    def validate_token_request(self) -> None: ...
    def create_token_response(self) -> _ServerResponse: ...

class AuthorizationEndpointMixin:
    RESPONSE_TYPES: Collection[str]
    ERROR_RESPONSE_FRAGMENT: bool
    @classmethod
    def check_authorization_endpoint(cls, request: OAuth2Request) -> bool: ...
    @staticmethod
    def validate_authorization_redirect_uri(request: OAuth2Request, client: ClientMixin) -> str: ...
    @staticmethod
    def validate_no_multiple_request_parameter(request: OAuth2Request): ...
    redirect_uri: str
    def validate_consent_request(self) -> str: ...
    def validate_authorization_request(self) -> str: ...
    def create_authorization_response(self, redirect_uri: str, grant_user) -> _ServerResponse: ...
