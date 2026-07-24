from _typeshed import Incomplete
from collections.abc import Mapping
from typing import overload
from typing_extensions import deprecated

from authlib.oauth2.rfc6749 import ClientMixin

class OAuth2Payload:
    @property
    def data(self) -> dict[str, str]: ...
    @property
    def datalist(self) -> dict[str, list[Incomplete]]: ...
    @property
    def client_id(self) -> str: ...
    @property
    def response_type(self) -> str: ...
    @property
    def grant_type(self) -> str: ...
    @property
    def redirect_uri(self) -> str: ...
    @property
    def scope(self) -> str: ...
    @property
    def state(self) -> str | None: ...

class BasicOAuth2Payload(OAuth2Payload):
    def __init__(self, payload: dict[str, str]) -> None: ...
    @property
    def data(self) -> dict[str, str]: ...
    @property
    def datalist(self) -> dict[str, list[Incomplete]]: ...

class OAuth2Request(OAuth2Payload):
    method: str
    uri: str
    headers: Mapping[str, str] | None
    payload: OAuth2Payload | None
    client: ClientMixin | None
    auth_method: str | None
    user: Incomplete | None
    authorization_code: Incomplete | None
    refresh_token: Incomplete | None
    credential: Incomplete | None
    @overload
    def __init__(self, method: str, uri: str, body: None = None, headers: Mapping[str, str] | None = None) -> None: ...
    @overload
    @deprecated("The `body` parameter in OAuth2Request is deprecated. Use the payload system instead.")
    def __init__(self, method: str, uri: str, body, headers: Mapping[str, str] | None = None) -> None: ...
    @property
    def args(self) -> dict[str, str | None]: ...
    @property
    def form(self) -> dict[str, str]: ...
    @property
    @deprecated("'request.data' is deprecated in favor of 'request.payload.data'")
    def data(self) -> dict[str, str]: ...
    @property
    @deprecated("'request.datalist' is deprecated in favor of 'request.payload.datalist'")
    def datalist(self) -> dict[str, list[Incomplete]]: ...
    @property
    @deprecated("'request.client_id' is deprecated in favor of 'request.payload.client_id'")
    def client_id(self) -> str: ...
    @property
    @deprecated("'request.response_type' is deprecated in favor of 'request.payload.response_type'")
    def response_type(self) -> str: ...
    @property
    @deprecated("'request.grant_type' is deprecated in favor of 'request.payload.grant_type'")
    def grant_type(self) -> str: ...
    @property
    @deprecated("'request.redirect_uri' is deprecated in favor of 'request.payload.redirect_uri'")
    def redirect_uri(self) -> str: ...
    @property
    @deprecated("'request.scope' is deprecated in favor of 'request.payload.scope'")
    def scope(self) -> str: ...
    @property
    @deprecated("'request.state' is deprecated in favor of 'request.payload.state'")
    def state(self) -> str | None: ...
    @property
    @deprecated("'request.body' is deprecated. Use the payload system instead.")
    def body(self): ...

class JsonPayload:
    @property
    def data(self): ...

class JsonRequest:
    method: str
    uri: str
    payload: JsonPayload | None
    headers: Mapping[str, str]
    def __init__(self, method: str, uri: str, headers: Mapping[str, str] | None = None) -> None: ...
    @property
    @deprecated("'request.data' is deprecated in favor of 'request.payload.data'")
    def data(self): ...
