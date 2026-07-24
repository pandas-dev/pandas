from _typeshed import Incomplete
from collections.abc import Iterable
from logging import Logger

from .base import BaseEndpoint

log: Logger

class MetadataEndpoint(BaseEndpoint):
    raise_errors: bool
    endpoints: Iterable[BaseEndpoint]
    initial_claims: dict[str, Incomplete]
    claims: dict[str, Incomplete]
    def __init__(
        self, endpoints: Iterable[BaseEndpoint], claims: dict[str, Incomplete] = {}, raise_errors: bool = True
    ) -> None: ...
    def create_metadata_response(
        self, uri: str, http_method: str = "GET", body: str | None = None, headers: dict[str, str] | None = None
    ) -> tuple[dict[str, str], str, int]: ...
    def validate_metadata(
        self, array, key, is_required: bool = False, is_list: bool = False, is_url: bool = False, is_issuer: bool = False
    ) -> None: ...
    def validate_metadata_token(self, claims, endpoint: BaseEndpoint) -> None: ...
    def validate_metadata_authorization(self, claims, endpoint: BaseEndpoint): ...
    def validate_metadata_revocation(self, claims, endpoint: BaseEndpoint) -> None: ...
    def validate_metadata_introspection(self, claims, endpoint: BaseEndpoint) -> None: ...
    def validate_metadata_server(self) -> dict[str, Incomplete]: ...
