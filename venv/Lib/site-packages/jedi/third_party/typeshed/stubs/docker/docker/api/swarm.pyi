import logging
from typing import Any, Literal, TypedDict, type_check_only
from typing_extensions import TypeAlias

from docker.types.services import DriverConfig
from docker.types.swarm import SwarmExternalCA, SwarmSpec

log: logging.Logger

@type_check_only
class _HasId(TypedDict):
    Id: str

@type_check_only
class _HasID(TypedDict):
    ID: str

_Node: TypeAlias = _HasId | _HasID | str

@type_check_only
class _NodeSpec(TypedDict, total=False):
    Name: str
    Labels: dict[str, str]
    Role: Literal["worker", "manager"]
    Availability: Literal["active", "pause", "drain"]

@type_check_only
class _UnlockKeyResponse(TypedDict):
    UnlockKey: str

class SwarmApiMixin:
    def create_swarm_spec(
        self,
        task_history_retention_limit: int | None = None,
        snapshot_interval: int | None = None,
        keep_old_snapshots: int | None = None,
        log_entries_for_slow_followers: int | None = None,
        heartbeat_tick: int | None = None,
        election_tick: int | None = None,
        dispatcher_heartbeat_period: int | None = None,
        node_cert_expiry: int | None = None,
        external_ca: SwarmExternalCA | None = None,
        external_cas: list[SwarmExternalCA] | None = None,
        name: str | None = None,
        labels: dict[str, str] | None = None,
        signing_ca_cert: str | None = None,
        signing_ca_key: str | None = None,
        ca_force_rotate: int | None = None,
        autolock_managers: bool | None = None,
        log_driver: DriverConfig | None = None,
    ) -> SwarmSpec: ...
    def get_unlock_key(self) -> _UnlockKeyResponse: ...
    def init_swarm(
        self,
        advertise_addr: str | None = None,
        listen_addr: str = "0.0.0.0:2377",
        force_new_cluster: bool = False,
        swarm_spec: dict[str, Any] | None = None,  # Any: arbitrary SwarmSpec configuration body
        default_addr_pool: list[str] | None = None,
        subnet_size: int | None = None,
        data_path_addr: str | None = None,
        data_path_port: int | None = None,
    ) -> str: ...
    def inspect_swarm(self) -> dict[str, Any]: ...  # Any: deeply nested ClusterInfo + JoinTokens
    def inspect_node(self, node_id: _Node) -> dict[str, Any]: ...  # Any: deeply nested Node object
    def join_swarm(
        self,
        remote_addrs: list[str],
        join_token: str,
        listen_addr: str = "0.0.0.0:2377",
        advertise_addr: str | None = None,
        data_path_addr: str | None = None,
    ) -> Literal[True]: ...
    def leave_swarm(self, force: bool = False) -> Literal[True]: ...
    def nodes(self, filters: dict[str, Any] | None = None) -> list[dict[str, Any]]: ...  # Any: filter values + Node response
    def remove_node(self, node_id: _Node, force: bool = False) -> Literal[True]: ...
    def unlock_swarm(self, key: str | _UnlockKeyResponse) -> Literal[True]: ...
    def update_node(self, node_id: _Node, version: int, node_spec: _NodeSpec | None = None) -> Literal[True]: ...
    def update_swarm(
        self,
        version: int,
        swarm_spec: dict[str, Any] | None = None,  # Any: arbitrary SwarmSpec configuration body
        rotate_worker_token: bool = False,
        rotate_manager_token: bool = False,
        rotate_manager_unlock_key: bool = False,
    ) -> Literal[True]: ...
