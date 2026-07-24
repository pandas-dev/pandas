from collections.abc import Iterable
from typing import Any, Literal, TypedDict, type_check_only
from typing_extensions import TypeAlias

from docker.types import IPAMConfig

@type_check_only
class _HasId(TypedDict):
    Id: str

@type_check_only
class _HasID(TypedDict):
    ID: str

_Network: TypeAlias = _HasId | _HasID | str
_Container: TypeAlias = _HasId | _HasID | str

class NetworkApiMixin:
    def networks(self, names=None, ids=None, filters=None): ...
    def create_network(
        self,
        name: str,
        driver: str | None = None,
        options: dict[str, Any] | None = None,
        ipam: IPAMConfig | None = None,
        check_duplicate: bool | None = None,
        internal: bool = False,
        labels: dict[str, Any] | None = None,
        enable_ipv6: bool = False,
        attachable: bool | None = None,
        scope: Literal["local", "global", "swarm"] | None = None,
        ingress: bool | None = None,
    ) -> dict[str, str]: ...
    def prune_networks(self, filters=None): ...
    def remove_network(self, net_id: _Network) -> None: ...
    def inspect_network(
        self, net_id: _Network, verbose: bool | None = None, scope: Literal["local", "global", "swarm"] | None = None
    ): ...
    def connect_container_to_network(
        self,
        container: _Container,
        net_id: str,
        ipv4_address=None,
        ipv6_address=None,
        aliases=None,
        links: dict[str, str] | dict[str, None] | dict[str, str | None] | Iterable[tuple[str, str | None]] | None = None,
        link_local_ips=None,
        driver_opt=None,
        mac_address=None,
    ) -> None: ...
    def disconnect_container_from_network(self, container: _Container, net_id: str, force: bool = False) -> None: ...
