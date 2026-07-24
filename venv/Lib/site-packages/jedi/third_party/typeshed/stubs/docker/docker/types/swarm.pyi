from typing import Any

from .services import DriverConfig

class SwarmSpec(dict[str, Any]):
    def __init__(
        self,
        version: str,
        task_history_retention_limit: int | None = None,
        snapshot_interval: int | None = None,
        keep_old_snapshots: int | None = None,
        log_entries_for_slow_followers: int | None = None,
        heartbeat_tick: int | None = None,
        election_tick: int | None = None,
        dispatcher_heartbeat_period: int | None = None,
        node_cert_expiry: int | None = None,
        external_cas: list[SwarmExternalCA] | None = None,
        name: str | None = None,
        labels: dict[str, str] | None = None,
        signing_ca_cert: str | None = None,
        signing_ca_key: str | None = None,
        ca_force_rotate: int | None = None,
        autolock_managers: bool | None = None,
        log_driver: DriverConfig | None = None,
    ) -> None: ...

class SwarmExternalCA(dict[str, Any]):
    def __init__(
        self, url: str, protocol: str | None = None, options: dict[str, str] | None = None, ca_cert: str | None = None
    ) -> None: ...
