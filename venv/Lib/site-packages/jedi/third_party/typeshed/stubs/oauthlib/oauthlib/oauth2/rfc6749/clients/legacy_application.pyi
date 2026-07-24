from _typeshed import Incomplete
from collections.abc import Callable

from .base import Client, _TokenPlacement

class LegacyApplicationClient(Client):
    grant_type: str
    def __init__(
        self,
        client_id: str,
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
    def prepare_request_body(
        self,
        username: str,
        password: str,
        body: str = "",
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        include_client_id: bool = False,
        *,
        code_verifier: str | None = None,
        client_id: str | None = None,
        client_secret: str | None = None,
        code: str | None = None,
        redirect_uri: str | None = None,
        **kwargs,
    ) -> str: ...
