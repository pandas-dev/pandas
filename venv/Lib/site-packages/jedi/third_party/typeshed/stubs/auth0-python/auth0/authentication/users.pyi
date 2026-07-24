from _typeshed import Incomplete

from auth0.rest import RestClient
from auth0.types import TimeoutType

class Users:
    domain: str
    protocol: str
    client: RestClient
    def __init__(self, domain: str, telemetry: bool = True, timeout: TimeoutType = 5.0, protocol: str = "https") -> None: ...
    def userinfo(self, access_token: str) -> dict[str, Incomplete]: ...
