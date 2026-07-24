from _typeshed import Incomplete
from collections.abc import Callable

from .base import Client, _TokenPlacement

class WebApplicationClient(Client):
    grant_type: str
    code: str | None
    def __init__(
        self,
        client_id: str,
        code: str | None = None,
        *,
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
    def prepare_request_uri(
        self,
        uri: str,
        redirect_uri: str | None = None,
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        state: str | None = None,
        code_challenge: str | None = None,
        code_challenge_method: str | None = "plain",
        **kwargs,
    ) -> str: ...
    def prepare_request_body(
        self,
        code: str | None = None,
        redirect_uri: str | None = None,
        body: str = "",
        include_client_id: bool = True,
        code_verifier: str | None = None,
        *,
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        client_id: str | None = None,
        client_secret: str | None = None,
        **kwargs,
    ) -> str: ...
    def parse_request_uri_response(self, uri: str, state: str | None = None) -> dict[str, str]: ...
