from collections.abc import Iterable
from typing import Any

from .resource import Model

class Swarm(Model):
    id_attribute: str
    def __init__(self, *args, **kwargs) -> None: ...
    @property
    def version(self) -> str | None: ...
    def get_unlock_key(self) -> dict[str, Any]: ...
    def init(
        self,
        advertise_addr: str | None = None,
        listen_addr: str = "0.0.0.0:2377",
        force_new_cluster: bool = False,
        default_addr_pool: Iterable[str] | None = None,
        subnet_size: int | None = None,
        data_path_addr: str | None = None,
        data_path_port: int | None = None,
        **kwargs,
    ) -> str: ...
    def join(self, *args, **kwargs) -> bool: ...
    def leave(self, *args, **kwargs) -> bool: ...
    def reload(self) -> None: ...
    def unlock(self, key: str) -> bool: ...
    def update(
        self,
        rotate_worker_token: bool = False,
        rotate_manager_token: bool = False,
        rotate_manager_unlock_key: bool = False,
        **kwargs,
    ) -> bool: ...
