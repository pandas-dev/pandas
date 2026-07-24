from collections.abc import Iterable, Mapping
from typing import Any, Final, Literal

from docker._types import ContainerWeightDevice

from .. import errors
from .base import DictType
from .healthcheck import Healthcheck
from .networks import NetworkingConfig
from .services import Mount

class LogConfigTypesEnum:
    JSON: Final = "json-file"
    SYSLOG: Final = "syslog"
    JOURNALD: Final = "journald"
    GELF: Final = "gelf"
    FLUENTD: Final = "fluentd"
    NONE: Final = "none"

class LogConfig(DictType[Any]):
    types: type[LogConfigTypesEnum]
    def __init__(
        self, *, type: str = ..., Type: str = ..., config: dict[str, str] = ..., Config: dict[str, str] = ...
    ) -> None: ...
    @property
    def type(self) -> str: ...
    @type.setter
    def type(self, value: str) -> None: ...
    @property
    def config(self) -> dict[str, str]: ...
    def set_config_value(self, key: str, value: str) -> None: ...
    def unset_config(self, key: str) -> None: ...

class Ulimit(DictType[Any]):
    def __init__(
        self, *, name: str = ..., Name: str = ..., soft: int = ..., Soft: int = ..., hard: int = ..., Hard: int = ...
    ) -> None: ...
    @property
    def name(self) -> str: ...
    @name.setter
    def name(self, value: str) -> None: ...
    @property
    def soft(self) -> int | None: ...
    @soft.setter
    def soft(self, value: int | None) -> None: ...
    @property
    def hard(self) -> int | None: ...
    @hard.setter
    def hard(self, value: int | None) -> None: ...

class DeviceRequest(DictType[Any]):
    def __init__(
        self,
        *,
        driver: str = ...,
        Driver: str = ...,
        count: int = ...,
        Count: int = ...,
        device_ids: list[str] = ...,
        DeviceIDs: list[str] = ...,
        capabilities: list[list[str]] = ...,
        Capabilities: list[list[str]] = ...,
        options: dict[str, str] = ...,
        Options: dict[str, str] = ...,
    ) -> None: ...
    @property
    def driver(self) -> str: ...
    @driver.setter
    def driver(self, value: str) -> None: ...
    @property
    def count(self) -> int: ...
    @count.setter
    def count(self, value: int) -> None: ...
    @property
    def device_ids(self) -> list[str]: ...
    @device_ids.setter
    def device_ids(self, value: list[str]) -> None: ...
    @property
    def capabilities(self) -> list[list[str]]: ...
    @capabilities.setter
    def capabilities(self, value: list[list[str]]) -> None: ...
    @property
    def options(self) -> dict[str, str]: ...
    @options.setter
    def options(self, value: dict[str, str]) -> None: ...

class HostConfig(dict[str, Any]):
    def __init__(
        self,
        version: str,
        binds: dict[str, Mapping[str, str]] | list[str] | None = None,
        port_bindings: Mapping[int | str, Any] | None = None,  # Any: int, str, tuple, dict, or list
        lxc_conf: dict[str, str] | list[dict[str, str]] | None = None,
        publish_all_ports: bool = False,
        links: dict[str, str] | dict[str, None] | dict[str, str | None] | Iterable[tuple[str, str | None]] | None = None,
        privileged: bool = False,
        dns: list[str] | None = None,
        dns_search: list[str] | None = None,
        volumes_from: list[str] | None = None,
        network_mode: str | None = None,
        restart_policy: Mapping[str, str | int] | None = None,
        cap_add: list[str] | None = None,
        cap_drop: list[str] | None = None,
        devices: list[str] | None = None,
        extra_hosts: dict[str, str] | list[str] | None = None,
        read_only: bool | None = None,
        pid_mode: str | None = None,
        ipc_mode: str | None = None,
        security_opt: list[str] | None = None,
        ulimits: list[Ulimit] | None = None,
        log_config: LogConfig | None = None,
        mem_limit: str | int | None = None,
        memswap_limit: str | int | None = None,
        mem_reservation: str | int | None = None,
        kernel_memory: str | int | None = None,
        mem_swappiness: int | None = None,
        cgroup_parent: str | None = None,
        group_add: Iterable[str | int] | None = None,
        cpu_quota: int | None = None,
        cpu_period: int | None = None,
        blkio_weight: int | None = None,
        blkio_weight_device: list[ContainerWeightDevice] | None = None,
        device_read_bps: list[Mapping[str, str | int]] | None = None,
        device_write_bps: list[Mapping[str, str | int]] | None = None,
        device_read_iops: list[Mapping[str, str | int]] | None = None,
        device_write_iops: list[Mapping[str, str | int]] | None = None,
        oom_kill_disable: bool = False,
        shm_size: str | int | None = None,
        sysctls: dict[str, str] | None = None,
        tmpfs: dict[str, str] | None = None,
        oom_score_adj: int | None = None,
        dns_opt: list[str] | None = None,
        cpu_shares: int | None = None,
        cpuset_cpus: str | None = None,
        userns_mode: str | None = None,
        uts_mode: str | None = None,
        pids_limit: int | None = None,
        isolation: str | None = None,
        auto_remove: bool = False,
        storage_opt: dict[str, str] | None = None,
        init: bool | None = None,
        init_path: str | None = None,
        volume_driver: str | None = None,
        cpu_count: int | None = None,
        cpu_percent: int | None = None,
        nano_cpus: int | None = None,
        cpuset_mems: str | None = None,
        runtime: str | None = None,
        mounts: list[Mount] | None = None,
        cpu_rt_period: int | None = None,
        cpu_rt_runtime: int | None = None,
        device_cgroup_rules: list[str] | None = None,
        device_requests: list[DeviceRequest] | None = None,
        cgroupns: Literal["private", "host"] | None = None,
    ) -> None: ...

def host_config_type_error(param: str, param_value: object, expected: str) -> TypeError: ...
def host_config_version_error(param: str, version: str, less_than: bool = True) -> errors.InvalidVersion: ...
def host_config_value_error(param: str, param_value: object) -> ValueError: ...
def host_config_incompatible_error(param: str, param_value: str, incompatible_param: str) -> errors.InvalidArgument: ...

class ContainerConfig(dict[str, Any]):
    def __init__(
        self,
        version: str,
        image: str,
        command: str | list[str],
        hostname: str | None = None,
        user: str | int | None = None,
        detach: bool = False,
        stdin_open: bool = False,
        tty: bool = False,
        # list is invariant, enumerating all possible union combination would be too complex for:
        # list[str | int | tuple[int | str, str] | tuple[int | str, ...]]
        ports: dict[str, dict[str, str]] | list[Any] | None = None,
        environment: dict[str, str] | list[str] | None = None,
        volumes: str | list[str] | None = None,
        network_disabled: bool = False,
        entrypoint: str | list[str] | None = None,
        working_dir: str | None = None,
        domainname: str | None = None,
        host_config: HostConfig | None = None,
        mac_address: str | None = None,
        labels: dict[str, str] | list[str] | None = None,
        stop_signal: str | None = None,
        networking_config: NetworkingConfig | None = None,
        healthcheck: Healthcheck | None = None,
        stop_timeout: int | None = None,
        runtime: str | None = None,
    ) -> None: ...
