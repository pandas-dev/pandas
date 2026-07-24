from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class Connections:
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
        strategy: str | None = None,
        fields: list[str] | None = None,
        include_fields: bool = True,
        page: int | None = None,
        per_page: int | None = None,
        extra_params: dict[str, Incomplete] | None = None,
        name: str | None = None,
    ) -> list[dict[str, Incomplete]]: ...
    async def all_async(
        self,
        strategy: str | None = None,
        fields: list[str] | None = None,
        include_fields: bool = True,
        page: int | None = None,
        per_page: int | None = None,
        extra_params: dict[str, Incomplete] | None = None,
        name: str | None = None,
    ) -> list[dict[str, Incomplete]]: ...
    def get(self, id: str, fields: list[str] | None = None, include_fields: bool = True) -> dict[str, Incomplete]: ...
    async def get_async(self, id: str, fields: list[str] | None = None, include_fields: bool = True) -> dict[str, Incomplete]: ...
    def delete(self, id: str): ...
    async def delete_async(self, id: str): ...
    def update(self, id: str, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def update_async(self, id: str, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    def create(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def create_async(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    def delete_user_by_email(self, id: str, email: str): ...
    async def delete_user_by_email_async(self, id: str, email: str): ...
