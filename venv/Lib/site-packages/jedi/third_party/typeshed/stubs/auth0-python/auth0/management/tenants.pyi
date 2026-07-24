from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class Tenants:
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
    def get(self, fields: list[str] | None = None, include_fields: bool = True) -> dict[str, Incomplete]: ...
    async def get_async(self, fields: list[str] | None = None, include_fields: bool = True) -> dict[str, Incomplete]: ...
    def update(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def update_async(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
