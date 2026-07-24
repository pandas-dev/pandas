from _typeshed import Incomplete
from collections.abc import Callable
from typing import Literal

from .tokens import OAuth2Token

def prepare_grant_uri(
    uri: str,
    client_id: str,
    response_type: Literal["code", "token"],
    redirect_uri: str | None = None,
    scope: str | set[object] | tuple[object] | list[object] | None = None,
    state: str | None = None,
    code_challenge: str | None = None,
    code_challenge_method: str | None = "plain",
    **kwargs,
) -> str: ...
def prepare_token_request(
    grant_type: str,
    body: str = "",
    include_client_id: bool = True,
    code_verifier: str | None = None,
    *,
    scope: str | set[object] | tuple[object] | list[object] | None = None,
    client_id: str | None = None,
    client_secret: str | None = None,
    code: str | None = None,
    redirect_uri: str | None = None,
    **kwargs,
) -> str: ...
def prepare_token_revocation_request(
    url: str,
    token: str,
    token_type_hint: Literal["access_token", "refresh_token"] | None = "access_token",
    callback: Callable[[Incomplete], Incomplete] | None = None,
    body: str = "",
    **kwargs,
) -> tuple[str, dict[str, str], str]: ...
def parse_authorization_code_response(uri: str, state: str | None = None) -> dict[str, str]: ...
def parse_implicit_response(
    uri: str, state: str | None = None, scope: str | set[object] | tuple[object] | list[object] | None = None
) -> OAuth2Token: ...
def parse_token_response(
    body: str | bytes | bytearray, scope: str | set[object] | tuple[object] | list[object] | None = None
) -> OAuth2Token: ...
def validate_token_parameters(params: dict[str, Incomplete]) -> None: ...
def parse_expires(params: dict[str, Incomplete]) -> tuple[int | None, float | None, float | None]: ...
