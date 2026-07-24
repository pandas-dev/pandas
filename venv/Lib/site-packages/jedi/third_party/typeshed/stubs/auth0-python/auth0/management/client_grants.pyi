from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class ClientGrants:
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
        audience: str | None = None,
        page: int | None = None,
        per_page: int | None = None,
        include_totals: bool = False,
        client_id: str | None = None,
        allow_any_organization: bool | None = None,
    ) -> dict[str, Incomplete]: ...
    async def all_async(
        self,
        audience: str | None = None,
        page: int | None = None,
        per_page: int | None = None,
        include_totals: bool = False,
        client_id: str | None = None,
        allow_any_organization: bool | None = None,
    ) -> dict[str, Incomplete]: ...
    def create(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def create_async(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    def delete(self, id: str): ...
    async def delete_async(self, id: str): ...
    def update(self, id: str, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def update_async(self, id: str, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    def get_organizations(
        self,
        id: str,
        page: int | None = None,
        per_page: int | None = None,
        include_totals: bool = False,
        from_param: str | None = None,
        take: int | None = None,
    ) -> dict[str, Incomplete]: ...
    async def get_organizations_async(
        self,
        id: str,
        page: int | None = None,
        per_page: int | None = None,
        include_totals: bool = False,
        from_param: str | None = None,
        take: int | None = None,
    ) -> dict[str, Incomplete]: ...
