from _typeshed import Incomplete
from logging import Logger
from typing import Any, Literal, Protocol, TypedDict, overload, type_check_only
from typing_extensions import TypeAlias

import requests
from oauthlib.oauth2 import Client
from requests.cookies import RequestsCookieJar

_Token: TypeAlias = dict[str, Incomplete]  # oauthlib.oauth2.Client.token

@type_check_only
class _AccessTokenResponseHook(Protocol):
    def __call__(self, response: requests.Response, /) -> requests.Response: ...

@type_check_only
class _RefreshTokenResponseHook(Protocol):
    def __call__(self, response: requests.Response, /) -> requests.Response: ...

@type_check_only
class _ProtectedRequestHook(Protocol):
    def __call__(self, url, headers, data, /) -> tuple[Incomplete, Incomplete, Incomplete]: ...

@type_check_only
class _ComplianceHooks(TypedDict):
    access_token_response: set[_AccessTokenResponseHook]
    refresh_token_response: set[_RefreshTokenResponseHook]
    protected_request: set[_ProtectedRequestHook]

log: Logger

class TokenUpdated(Warning):
    token: Incomplete
    def __init__(self, token) -> None: ...

class OAuth2Session(requests.Session):
    redirect_uri: Incomplete
    state: Incomplete
    auto_refresh_url: str | None
    auto_refresh_kwargs: dict[str, Any]
    token_updater: Incomplete
    compliance_hook: _ComplianceHooks
    def __init__(
        self,
        client_id=None,
        client: Client | None = None,
        auto_refresh_url: str | None = None,
        auto_refresh_kwargs: dict[str, Any] | None = None,
        scope=None,
        redirect_uri=None,
        token=None,
        state=None,
        token_updater=None,
        pkce=None,
        **kwargs,
    ) -> None: ...
    @property
    def scope(self) -> Incomplete | None: ...  # oauthlib.oauth2.Client.scope
    @scope.setter
    def scope(self, value: Incomplete | None) -> None: ...
    def new_state(self): ...
    @property
    def client_id(self) -> Incomplete | None: ...  # oauthlib.oauth2.Client.client_id
    @client_id.setter
    def client_id(self, value: Incomplete | None) -> None: ...
    @client_id.deleter
    def client_id(self) -> None: ...
    @property
    def token(self): ...  # oauthlib.oauth2.Client.token
    @token.setter
    def token(self, value) -> None: ...
    @property
    def access_token(self): ...  # oauthlib.oauth2.Client.access_token
    @access_token.setter
    def access_token(self, value) -> None: ...
    @access_token.deleter
    def access_token(self) -> None: ...
    @property
    def authorized(self) -> bool: ...
    def authorization_url(self, url: str, state=None, **kwargs) -> tuple[str, str]: ...
    def fetch_token(
        self,
        token_url: str,
        code=None,
        authorization_response=None,
        body: str = "",
        auth=None,
        username=None,
        password=None,
        method: str = "POST",
        force_querystring: bool = False,
        timeout=None,
        headers=None,
        verify: bool | None = None,
        proxies=None,
        include_client_id=None,
        client_secret=None,
        cert=None,
        **kwargs,
    ) -> _Token: ...
    def token_from_fragment(self, authorization_response: str) -> _Token: ...
    def refresh_token(
        self,
        token_url: str,
        refresh_token=None,
        body: str = "",
        auth=None,
        timeout=None,
        headers=None,
        verify: bool | None = None,
        proxies=None,
        **kwargs,
    ) -> _Token: ...
    def request(  # type: ignore[override]
        self,
        method: str | bytes,
        url: str | bytes,
        data: requests.sessions._Data | None = None,
        headers: requests.sessions._HeadersUpdateMapping | None = None,
        withhold_token: bool = False,
        client_id=None,
        client_secret=None,
        files: requests.sessions._Files | None = None,
        *,
        params: requests.sessions._Params | None = None,
        cookies: None | RequestsCookieJar | requests.sessions._TextMapping = None,
        auth: requests.sessions._Auth | None = None,
        timeout: requests.sessions._Timeout | None = None,
        allow_redirects: bool = True,
        proxies: requests.sessions._TextMapping | None = None,
        hooks: requests.sessions._HooksInput | None = None,
        stream: bool | None = None,
        verify: requests.sessions._Verify | None = None,
        cert: requests.sessions._Cert | None = None,
        json=None,
    ) -> requests.Response: ...
    @overload
    def register_compliance_hook(self, hook_type: Literal["access_token_response"], hook: _AccessTokenResponseHook) -> None: ...
    @overload
    def register_compliance_hook(self, hook_type: Literal["refresh_token_response"], hook: _RefreshTokenResponseHook) -> None: ...
    @overload
    def register_compliance_hook(self, hook_type: Literal["protected_request"], hook: _ProtectedRequestHook) -> None: ...
