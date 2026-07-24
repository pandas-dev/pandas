from _typeshed import ConvertibleToInt, Incomplete
from collections.abc import Callable
from typing import Final, Literal
from typing_extensions import TypeAlias

from oauthlib.common import _HTTPMethod
from oauthlib.oauth2.rfc6749.tokens import OAuth2Token

_TokenPlacement: TypeAlias = Literal["auth_header", "query", "body"]

AUTH_HEADER: Final[_TokenPlacement]
URI_QUERY: Final[_TokenPlacement]
BODY: Final[_TokenPlacement]
FORM_ENC_HEADERS: Final[dict[str, str]]

class Client:
    refresh_token_key: str
    client_id: str
    default_token_placement: _TokenPlacement
    token_type: str
    access_token: str | None
    refresh_token: str | None
    mac_key: str | bytes | bytearray | None
    mac_algorithm: str | None
    token: dict[str, Incomplete]
    scope: str | set[object] | tuple[object] | list[object]
    state_generator: Callable[[], str]
    state: str | None
    redirect_url: str | None
    code: str | None
    expires_in: ConvertibleToInt | None
    code_verifier: str | None
    code_challenge: str | None
    code_challenge_method: str | None
    def __init__(
        self,
        client_id: str,
        default_token_placement: _TokenPlacement = "auth_header",
        token_type: str = "Bearer",
        access_token: str | None = None,
        refresh_token: str | None = None,
        mac_key: str | bytes | bytearray | None = None,
        mac_algorithm: str | None = None,
        token: dict[str, Incomplete] | None = None,
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        state: str | None = None,
        redirect_url: str | None = None,
        state_generator: Callable[[], str] = ...,
        code_verifier: str | None = None,
        code_challenge: str | None = None,
        code_challenge_method: str | None = None,
        **kwargs,
    ) -> None: ...
    @property
    def token_types(
        self,
    ) -> dict[
        Literal["Bearer", "MAC"],
        Callable[
            [str, str, str | None, dict[str, str] | None, str | None, Incomplete], tuple[str, dict[str, str] | None, str | None]
        ],
    ]: ...
    def prepare_request_uri(self, *args, **kwargs) -> str: ...
    def prepare_request_body(self, *args, **kwargs) -> str: ...
    def parse_request_uri_response(self, *args, **kwargs) -> dict[str, str]: ...
    def add_token(
        self,
        uri: str,
        http_method: _HTTPMethod = "GET",
        body: str | None = None,
        headers: dict[str, str] | None = None,
        token_placement: _TokenPlacement | None = None,
        **kwargs,
    ) -> tuple[str, dict[str, str] | None, str | None]: ...
    def prepare_authorization_request(
        self,
        authorization_url: str,
        state: str | None = None,
        redirect_url: str | None = None,
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        **kwargs,
    ) -> tuple[str, dict[str, str], str]: ...
    def prepare_token_request(
        self,
        token_url: str,
        authorization_response: str | None = None,
        redirect_url: str | None = None,
        state: str | None = None,
        body: str = "",
        **kwargs,
    ) -> tuple[str, dict[str, str], str]: ...
    def prepare_refresh_token_request(
        self,
        token_url: str,
        refresh_token: str | None = None,
        body: str = "",
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        **kwargs,
    ) -> tuple[str, dict[str, str], str]: ...
    def prepare_token_revocation_request(
        self,
        revocation_url: str,
        token: str,
        token_type_hint: Literal["access_token", "refresh_token"] | None = "access_token",
        body: str = "",
        callback: Callable[[Incomplete], Incomplete] | None = None,
        **kwargs,
    ): ...
    def parse_request_body_response(
        self, body: str, scope: str | set[object] | tuple[object] | list[object] | None = None, **kwargs
    ) -> OAuth2Token: ...
    def prepare_refresh_body(
        self,
        body: str = "",
        refresh_token: str | None = None,
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        **kwargs,
    ) -> str: ...
    def create_code_verifier(self, length: int) -> str: ...
    def create_code_challenge(self, code_verifier: str, code_challenge_method: str | None = None) -> str: ...
    def populate_code_attributes(self, response: dict[str, Incomplete]) -> None: ...
    def populate_token_attributes(self, response: dict[str, Incomplete]) -> None: ...
