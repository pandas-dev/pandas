from collections.abc import Iterable, Mapping
from typing import Any, Final, Literal, TypedDict, TypeVar, overload, type_check_only

from .healthcheck import Healthcheck

_T = TypeVar("_T")

class TaskTemplate(dict[str, Any]):
    def __init__(
        self,
        container_spec: ContainerSpec,
        resources: Resources | None = None,
        restart_policy: RestartPolicy | None = None,
        placement: Placement | list[str] | None = None,
        log_driver: DriverConfig | None = None,
        networks: Iterable[str | NetworkAttachmentConfig] | None = None,
        force_update: int | None = None,
    ) -> None: ...
    @property
    def container_spec(self) -> ContainerSpec: ...
    @property
    def resources(self) -> Resources: ...
    @property
    def restart_policy(self) -> RestartPolicy: ...
    @property
    def placement(self) -> Placement: ...

class ContainerSpec(dict[str, Any]):
    def __init__(
        self,
        image: str,
        command: str | list[str] | None = None,
        args: list[str] | None = None,
        hostname: str | None = None,
        env: dict[str, str | bytes | None] | list[str] | None = None,
        workdir: str | None = None,
        user: str | None = None,
        labels: dict[str, str] | None = None,
        mounts: Iterable[str | Mount] | None = None,
        stop_grace_period: int | None = None,
        secrets: list[SecretReference] | None = None,
        tty: bool | None = None,
        groups: list[str] | None = None,
        open_stdin: bool | None = None,
        read_only: bool | None = None,
        stop_signal: str | None = None,
        healthcheck: Healthcheck | None = None,
        hosts: Mapping[str, str] | None = None,
        dns_config: DNSConfig | None = None,
        configs: list[ConfigReference] | None = None,
        privileges: Privileges | None = None,
        isolation: str | None = None,
        init: bool | None = None,
        cap_add: list[str] | None = None,
        cap_drop: list[str] | None = None,
        sysctls: dict[str, str] | None = None,
    ) -> None: ...

class Mount(dict[str, Any]):
    def __init__(
        self,
        target: str,
        source: str | None,
        type: Literal["bind", "volume", "tmpfs", "npipe"] = "volume",
        read_only: bool = False,
        consistency: Literal["default", "consistent", "cached", "delegated"] | None = None,
        propagation: str | None = None,
        no_copy: bool = False,
        labels: dict[str, str] | None = None,
        driver_config: DriverConfig | None = None,
        tmpfs_size: int | str | None = None,
        tmpfs_mode: int | None = None,
    ) -> None: ...
    @classmethod
    def parse_mount_string(cls, string: str) -> Mount: ...

@type_check_only
class _ResourceDict(TypedDict):
    Kind: str
    Value: int

class Resources(dict[str, Any]):
    def __init__(
        self,
        cpu_limit: int | None = None,
        mem_limit: int | None = None,
        cpu_reservation: int | None = None,
        mem_reservation: int | None = None,
        generic_resources: (
            dict[str, int | str] | list[dict[Literal["DiscreteResourceSpec", "NamedResourceSpec"], _ResourceDict]] | None
        ) = None,
    ) -> None: ...

class UpdateConfig(dict[str, Any]):
    def __init__(
        self,
        parallelism: int = 0,
        delay: int | None = None,
        failure_action: Literal["pause", "continue", "rollback"] = "continue",
        monitor: int | None = None,
        max_failure_ratio: float | None = None,
        order: Literal["start-first", "stop-first"] | None = None,
    ) -> None: ...

class RollbackConfig(UpdateConfig): ...

class RestartConditionTypesEnum:
    NONE: Final = "none"
    ON_FAILURE: Final = "on-failure"
    ANY: Final = "any"

class RestartPolicy(dict[str, Any]):
    condition_types: type[RestartConditionTypesEnum]
    def __init__(
        self, condition: Literal["none", "on-failure", "any"] = "none", delay: int = 0, max_attempts: int = 0, window: int = 0
    ) -> None: ...

class DriverConfig(dict[str, Any]):
    def __init__(self, name: str, options: dict[str, str] | None = None) -> None: ...

class EndpointSpec(dict[str, Any]):
    def __init__(
        self, mode: str | None = None, ports: Mapping[str, str | tuple[str | None, ...]] | list[dict[str, str]] | None = None
    ) -> None: ...

@overload
def convert_service_ports(ports: list[_T]) -> list[_T]: ...
@overload
def convert_service_ports(ports: Mapping[str, str | tuple[str | None, ...]]) -> list[dict[str, str]]: ...

class ServiceMode(dict[str, Any]):
    mode: Literal["replicated", "global", "ReplicatedJob", "GlobalJob"]
    def __init__(
        self,
        mode: Literal["replicated", "global", "replicated-job", "global-job"],
        replicas: int | None = None,
        concurrency: int | None = None,
    ) -> None: ...
    @property
    def replicas(self) -> int | None: ...

class SecretReference(dict[str, Any]):
    def __init__(
        self,
        secret_id: str,
        secret_name: str,
        filename: str | None = None,
        uid: str | None = None,
        gid: str | None = None,
        mode: int = 292,
    ) -> None: ...

class ConfigReference(dict[str, Any]):
    def __init__(
        self,
        config_id: str,
        config_name: str,
        filename: str | None = None,
        uid: str | None = None,
        gid: str | None = None,
        mode: int = 292,
    ) -> None: ...

class Placement(dict[str, Any]):
    def __init__(
        self,
        constraints: list[str] | None = None,
        preferences: Iterable[tuple[str, str] | PlacementPreference] | None = None,
        platforms: Iterable[tuple[str, str]] | None = None,
        maxreplicas: int | None = None,
    ) -> None: ...

class PlacementPreference(dict[str, Any]):
    def __init__(self, strategy: Literal["spread"], descriptor: str) -> None: ...

class DNSConfig(dict[str, Any]):
    def __init__(
        self, nameservers: list[str] | None = None, search: list[str] | None = None, options: list[str] | None = None
    ) -> None: ...

class Privileges(dict[str, Any]):
    def __init__(
        self,
        credentialspec_file: str | None = None,
        credentialspec_registry: str | None = None,
        selinux_disable: bool | None = None,
        selinux_user: str | None = None,
        selinux_role: str | None = None,
        selinux_type: str | None = None,
        selinux_level: str | None = None,
    ) -> None: ...

class NetworkAttachmentConfig(dict[str, Any]):
    def __init__(self, target: str, aliases: list[str] | None = None, options: dict[str, str] | None = None) -> None: ...
