from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class Tickets:
    domain: Incomplete
    protocol: Incomplete
    client: Incomplete
    def __init__(
        self,
        domain: str,
        token: str,
        telemetry: bool = True,
        timeout: TimeoutType = 5.0,
        protocol: str = "https",
        rest_options: RestClientOptions | None = None,
    ) -> None: ...
    def create_email_verification(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def create_email_verification_async(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    def create_pswd_change(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def create_pswd_change_async(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
