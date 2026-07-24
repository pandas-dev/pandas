from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class EmailTemplates:
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
    def create(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def create_async(self, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    def get(self, template_name: str) -> dict[str, Incomplete]: ...
    async def get_async(self, template_name: str) -> dict[str, Incomplete]: ...
    def update(self, template_name: str, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
    async def update_async(self, template_name: str, body: dict[str, Incomplete]) -> dict[str, Incomplete]: ...
