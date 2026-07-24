from _typeshed import Incomplete
from collections.abc import Callable

from oauthlib.oauth2.rfc6749.clients.base import Client, _TokenPlacement

class DeviceClient(Client):
    grant_type: str
    client_secret: str | None
    def __init__(
        self,
        client_id: str,
        *,
        client_secret: str | None = None,
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
        self, uri: str, scope: str | set[object] | tuple[object] | list[object] | None = None, **kwargs
    ) -> str: ...
    def prepare_request_body(
        self,
        device_code: str,
        body: str = "",
        scope: str | set[object] | tuple[object] | list[object] | None = None,
        include_client_id: bool = False,
        **kwargs,
    ) -> str: ...
