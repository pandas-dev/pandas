from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class Logs:
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
    def search(
        self,
        page: int = 0,
        per_page: int = 50,
        sort: str | None = None,
        q: str | None = None,
        include_totals: bool = True,
        fields: list[str] | None = None,
        from_param: str | None = None,
        take: int | None = None,
        include_fields: bool = True,
    ) -> dict[str, Incomplete]: ...
    async def search_async(
        self,
        page: int = 0,
        per_page: int = 50,
        sort: str | None = None,
        q: str | None = None,
        include_totals: bool = True,
        fields: list[str] | None = None,
        from_param: str | None = None,
        take: int | None = None,
        include_fields: bool = True,
    ) -> dict[str, Incomplete]: ...
    def get(self, id: str) -> dict[str, Incomplete]: ...
    async def get_async(self, id: str) -> dict[str, Incomplete]: ...
