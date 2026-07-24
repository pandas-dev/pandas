from _typeshed import Incomplete
from typing import Final

from auth0.rest import RestClient
from auth0.types import RequestData

UNKNOWN_ERROR: Final[str]

class AuthenticationBase:
    domain: str
    client_id: str
    client_secret: str | None
    client_assertion_signing_key: str | None
    client_assertion_signing_alg: str | None
    protocol: str
    client: RestClient
    def __init__(
        self,
        domain: str,
        client_id: str,
        client_secret: str | None = None,
        client_assertion_signing_key: str | None = None,
        client_assertion_signing_alg: str | None = None,
        telemetry: bool = True,
        timeout: float | tuple[float, float] = 5.0,
        protocol: str = "https",
    ) -> None: ...
    def post(self, url: str, data: RequestData | None = None, headers: dict[str, str] | None = None): ...
    def authenticated_post(self, url: str, data: dict[str, Incomplete], headers: dict[str, str] | None = None): ...
    def get(self, url: str, params: dict[str, Incomplete] | None = None, headers: dict[str, str] | None = None): ...
