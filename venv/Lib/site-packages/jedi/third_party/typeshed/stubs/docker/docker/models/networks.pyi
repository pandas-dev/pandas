from typing import Any, Literal

from docker.types import IPAMConfig

from .containers import Container
from .resource import Collection, Model

class Network(Model):
    @property
    def name(self) -> str | None: ...
    @property
    def containers(self) -> list[Container]: ...
    def connect(self, container: str | Container, *args, **kwargs) -> None: ...
    def disconnect(self, container: str | Container, force: bool = False) -> None: ...
    def remove(self) -> None: ...

class NetworkCollection(Collection[Network]):
    model: type[Network]
    def create(  # type: ignore[override]
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
    ) -> Network: ...
    def get(
        self, network_id: str, verbose: bool | None = None, scope: Literal["local", "global", "swarm"] | None = None
    ) -> Network: ...  # type: ignore[override]
    def list(self, *args, **kwargs) -> list[Network]: ...
    def prune(self, filters: dict[str, Any] | None = None) -> dict[str, Any]: ...
