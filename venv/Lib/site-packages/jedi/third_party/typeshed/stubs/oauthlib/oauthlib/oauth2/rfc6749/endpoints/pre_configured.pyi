from _typeshed import Unused
from collections.abc import Callable

from oauthlib.common import Request
from oauthlib.oauth2.rfc8628.grant_types import DeviceCodeGrant

from ..grant_types import (
    AuthorizationCodeGrant,
    ClientCredentialsGrant,
    ImplicitGrant,
    RefreshTokenGrant,
    ResourceOwnerPasswordCredentialsGrant,
)
from ..request_validator import RequestValidator
from ..tokens import BearerToken
from .authorization import AuthorizationEndpoint
from .introspect import IntrospectEndpoint
from .resource import ResourceEndpoint
from .revocation import RevocationEndpoint
from .token import TokenEndpoint

class Server(AuthorizationEndpoint, IntrospectEndpoint, TokenEndpoint, ResourceEndpoint, RevocationEndpoint):
    auth_grant: AuthorizationCodeGrant
    implicit_grant: ImplicitGrant
    password_grant: ResourceOwnerPasswordCredentialsGrant
    credentials_grant: ClientCredentialsGrant
    refresh_grant: RefreshTokenGrant
    device_code_grant: DeviceCodeGrant
    bearer: BearerToken
    def __init__(
        self,
        request_validator: RequestValidator,
        token_expires_in: int | Callable[[Request], int] | None = None,
        token_generator: Callable[[Request], str] | None = None,
        refresh_token_generator: Callable[[Request], str] | None = None,
        *args: Unused,
    ) -> None: ...

class WebApplicationServer(AuthorizationEndpoint, IntrospectEndpoint, TokenEndpoint, ResourceEndpoint, RevocationEndpoint):
    auth_grant: AuthorizationCodeGrant
    refresh_grant: RefreshTokenGrant
    bearer: BearerToken
    def __init__(
        self,
        request_validator: RequestValidator,
        token_generator: Callable[[Request], str] | None = None,
        token_expires_in: int | Callable[[Request], int] | None = None,
        refresh_token_generator: Callable[[Request], str] | None = None,
    ) -> None: ...

class MobileApplicationServer(AuthorizationEndpoint, IntrospectEndpoint, ResourceEndpoint, RevocationEndpoint):
    implicit_grant: ImplicitGrant
    bearer: BearerToken
    def __init__(
        self,
        request_validator: RequestValidator,
        token_generator: Callable[[Request], str] | None = None,
        token_expires_in: int | Callable[[Request], int] | None = None,
        refresh_token_generator: Callable[[Request], str] | None = None,
    ) -> None: ...

class LegacyApplicationServer(TokenEndpoint, IntrospectEndpoint, ResourceEndpoint, RevocationEndpoint):
    password_grant: ResourceOwnerPasswordCredentialsGrant
    refresh_grant: RefreshTokenGrant
    bearer: BearerToken
    def __init__(
        self,
        request_validator: RequestValidator,
        token_generator: Callable[[Request], str] | None = None,
        token_expires_in: int | Callable[[Request], int] | None = None,
        refresh_token_generator: Callable[[Request], str] | None = None,
    ) -> None: ...

class BackendApplicationServer(TokenEndpoint, IntrospectEndpoint, ResourceEndpoint, RevocationEndpoint):
    credentials_grant: ClientCredentialsGrant
    bearer: BearerToken
    def __init__(
        self,
        request_validator: RequestValidator,
        token_generator: Callable[[Request], str] | None = None,
        token_expires_in: int | Callable[[Request], int] | None = None,
        refresh_token_generator: Callable[[Request], str] | None = None,
    ) -> None: ...
