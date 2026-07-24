from _typeshed import Incomplete

from ..rest import RestClientOptions
from ..types import TimeoutType

class UsersByEmail:
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
    def search_users_by_email(
        self, email: str, fields: list[str] | None = None, include_fields: bool = True
    ) -> list[dict[str, Incomplete]]: ...
    async def search_users_by_email_async(
        self, email: str, fields: list[str] | None = None, include_fields: bool = True
    ) -> list[dict[str, Incomplete]]: ...
