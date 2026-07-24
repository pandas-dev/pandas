from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class DeviceCredentials:
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
    def get(
        self,
        user_id: str,
        client_id: str,
        type: str,
        fields: list[str] | None = None,
        include_fields: bool = True,
        page: int | None = None,
        per_page: int | None = None,
        include_totals: bool = False,
    ): ...
    async def get_async(
        self,
        user_id: str,
        client_id: str,
        type: str,
        fields: list[str] | None = None,
        include_fields: bool = True,
        page: int | None = None,
        per_page: int | None = None,
        include_totals: bool = False,
    ): ...
    def create(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def create_async(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    def delete(self, id: str): ...
    async def delete_async(self, id: str): ...
