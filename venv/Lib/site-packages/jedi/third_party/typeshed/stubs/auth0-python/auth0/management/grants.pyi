from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class Grants:
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
    def all(
        self,
        page: int | None = None,
        per_page: int | None = None,
        include_totals: bool = False,
        extra_params: dict[str, Incomplete] | None = None,
    ): ...
    async def all_async(
        self,
        page: int | None = None,
        per_page: int | None = None,
        include_totals: bool = False,
        extra_params: dict[str, Incomplete] | None = None,
    ): ...
    def delete(self, id: str): ...
    async def delete_async(self, id: str): ...
