import re
from _typeshed import Incomplete, SupportsLenAndGetItem
from collections.abc import Iterable, Mapping
from logging import Logger
from typing import Any, Final, Literal, TypeVar, overload
from typing_extensions import TypeAlias

_T = TypeVar("_T")
_V = TypeVar("_V")

_HTTPMethod: TypeAlias = Literal["CONNECT", "DELETE", "GET", "HEAD", "OPTIONS", "PATCH", "POST", "PUT", "TRACE"]

UNICODE_ASCII_CHARACTER_SET: Final[str]
CLIENT_ID_CHARACTER_SET: Final[str]
SANITIZE_PATTERN: Final[re.Pattern[str]]
INVALID_HEX_PATTERN: Final[re.Pattern[str]]
always_safe: Final[str]
log: Logger

def quote(s: str | bytes, safe: bytes = b"/") -> str: ...
def unquote(s: str | bytes) -> str: ...
def urlencode(params: Iterable[tuple[str | bytes, str | bytes]]) -> str: ...
def encode_params_utf8(params: Iterable[tuple[str | bytes, str | bytes]]) -> list[tuple[bytes, bytes]]: ...
def decode_params_utf8(params: Iterable[tuple[str | bytes, str | bytes]]) -> list[tuple[str, str]]: ...

urlencoded: Final[set[str]]

def urldecode(query: str | bytes) -> list[tuple[str, str]]: ...
def extract_params(raw: str | bytes | dict[str, str] | Iterable[tuple[str, str]]) -> list[tuple[str, str]] | None: ...
def generate_nonce() -> str: ...
def generate_timestamp() -> str: ...
def generate_token(length: int = 30, chars: SupportsLenAndGetItem[str] = ...) -> str: ...
def generate_signed_token(private_pem: str, request: Request) -> str: ...
def verify_signed_token(public_pem, token): ...
def generate_client_id(length: int = 30, chars: SupportsLenAndGetItem[str] = ...) -> str: ...
def add_params_to_qs(query: str, params: dict[str, str] | Iterable[tuple[str, str]]) -> str: ...
def add_params_to_uri(uri: str, params: dict[str, str] | Iterable[tuple[str, str]], fragment: bool = False) -> str: ...
def safe_string_equals(a: str, b: str) -> bool: ...
@overload
def to_unicode(data: str | bytes, encoding: str = "UTF-8") -> str: ...
@overload
def to_unicode(data: Mapping[str, _V] | Mapping[bytes, _V], encoding: str = "UTF-8") -> dict[str, _V]: ...
@overload
def to_unicode(data: _T, encoding: str = "UTF-8") -> _T: ...

class CaseInsensitiveDict(dict[str, Incomplete]):
    proxy: dict[str, str]
    def __init__(self, data: dict[str, Incomplete]) -> None: ...
    @overload
    def __contains__(self, k: str) -> bool: ...
    @overload
    def __contains__(self, k: object) -> bool: ...
    def __delitem__(self, k: str) -> None: ...
    def __getitem__(self, k: str): ...
    @overload
    def get(self, k: str, default: None = None) -> Incomplete | None: ...
    @overload
    def get(self, k: str, default): ...
    def __setitem__(self, k: str, v) -> None: ...
    def update(self, *args, **kwargs) -> None: ...

class Request:
    uri: str
    http_method: _HTTPMethod
    headers: CaseInsensitiveDict
    body: str | dict[str, str] | list[tuple[str, str]] | None
    decoded_body: list[tuple[str, str]] | None
    oauth_params: list[str]
    validator_log: dict[str, Any]  # value type depends on the key
    access_token: Incomplete | None
    client: Incomplete | None
    client_id: Incomplete | None
    client_secret: Incomplete | None
    code: Incomplete | None
    code_challenge: Incomplete | None
    code_challenge_method: Incomplete | None
    code_verifier: Incomplete | None
    extra_credentials: Incomplete | None
    grant_type: Incomplete | None
    redirect_uri: Incomplete | None
    refresh_token: Incomplete | None
    request_token: Incomplete | None
    response_type: Incomplete | None
    scope: Incomplete | None
    scopes: Incomplete | None
    state: Incomplete | None
    token: Incomplete | None
    user: Incomplete | None
    token_type_hint: Incomplete | None
    response_mode: Incomplete | None
    nonce: Incomplete | None
    display: Incomplete | None
    prompt: Incomplete | None
    claims: Incomplete | None
    max_age: Incomplete | None
    ui_locales: Incomplete | None
    id_token_hint: Incomplete | None
    login_hint: Incomplete | None
    acr_values: Incomplete | None
    def __init__(
        self,
        uri: str,
        http_method: _HTTPMethod = "GET",
        body: str | dict[str, str] | list[tuple[str, str]] | None = None,
        headers: Mapping[str, str] | None = None,
        encoding: str = "utf-8",
    ): ...
    def __getattr__(self, name: str) -> str | None: ...  # or raises AttributeError if attribute is not found
    @property
    def uri_query(self) -> str: ...
    @property
    def uri_query_params(self) -> list[tuple[str, str]]: ...
    @property
    def duplicate_params(self) -> list[str]: ...
