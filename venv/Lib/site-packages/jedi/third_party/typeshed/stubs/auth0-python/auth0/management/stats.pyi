from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class Stats:
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
    def active_users(self) -> int: ...
    async def active_users_async(self) -> int: ...
    def daily_stats(self, from_date: str | None = None, to_date: str | None = None) -> list[dict[str, Incomplete]]: ...
    async def daily_stats_async(
        self, from_date: str | None = None, to_date: str | None = None
    ) -> list[dict[str, Incomplete]]: ...
