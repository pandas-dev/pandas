from _typeshed import Incomplete
from collections.abc import Callable, Iterable
from itertools import chain
from logging import Logger
from typing import TypeVar
from typing_extensions import TypeAlias

from oauthlib.common import Request

from ..request_validator import RequestValidator
from ..tokens import TokenBase

log: Logger

_T = TypeVar("_T")
_AuthValidator: TypeAlias = Callable[[Request], dict[str, Incomplete]]
_TokenValidator: TypeAlias = Callable[[Request], None]
_CodeModifier: TypeAlias = Callable[[dict[str, str], TokenBase | None, Request | None], dict[str, str]]
_TokenModifier: TypeAlias = Callable[[dict[str, Incomplete], TokenBase | None, Request | None], dict[str, Incomplete]]

class ValidatorsContainer:
    pre_auth: Iterable[_AuthValidator]
    post_auth: Iterable[_AuthValidator]
    pre_token: Iterable[_TokenValidator]
    post_token: Iterable[_TokenValidator]
    def __init__(
        self,
        post_auth: Iterable[_AuthValidator],
        post_token: Iterable[_TokenValidator],
        pre_auth: Iterable[_AuthValidator],
        pre_token: Iterable[_TokenValidator],
    ) -> None: ...
    @property
    def all_pre(self) -> chain[_AuthValidator | _TokenValidator]: ...
    @property
    def all_post(self) -> chain[_AuthValidator | _TokenValidator]: ...

class GrantTypeBase:
    error_uri: str | None
    request_validator: RequestValidator | None
    default_response_mode: str
    refresh_token: bool
    response_types: list[str]
    def __init__(
        self,
        request_validator: RequestValidator | None = None,
        *,
        post_auth: Iterable[_AuthValidator] | None = None,
        post_token: Iterable[_TokenValidator] | None = None,
        pre_auth: Iterable[_AuthValidator] | None = None,
        pre_token: Iterable[_TokenValidator] | None = None,
        **kwargs,
    ) -> None: ...
    def register_response_type(self, response_type: str) -> None: ...
    def register_code_modifier(self, modifier: _CodeModifier) -> None: ...
    def register_token_modifier(self, modifier: _TokenModifier) -> None: ...
    def create_authorization_response(
        self, request: Request, token_handler: TokenBase
    ) -> tuple[dict[str, str], str | None, int | None]: ...
    def create_token_response(
        self, request: Request, token_handler: TokenBase
    ) -> tuple[dict[str, str], str | None, int | None]: ...
    def add_token(self, token: dict[str, _T], token_handler: TokenBase, request: Request) -> dict[str, _T]: ...
    def validate_grant_type(self, request: Request) -> None: ...
    def validate_scopes(self, request: Request) -> None: ...
    def prepare_authorization_response(
        self, request: Request, token: dict[str, Incomplete], headers: dict[str, str], body: str | None, status: int | None
    ) -> tuple[dict[str, str], str | None, int | None]: ...
