import datetime
from _io import _BufferedReaderStream
from collections.abc import Iterable, Iterator, Mapping
from socket import SocketIO
from typing import Any, Literal, NamedTuple, TypedDict, overload, type_check_only
from typing_extensions import NotRequired, override

from docker._types import ContainerWeightDevice, WaitContainerResponse
from docker.transport.sshconn import SSHSocket
from docker.types import EndpointConfig
from docker.types.containers import DeviceRequest, LogConfig, Ulimit
from docker.types.daemon import CancellableStream
from docker.types.services import Mount

from .images import Image
from .resource import Collection, Model

@type_check_only
class _RestartPolicy(TypedDict):
    MaximumRetryCount: NotRequired[int]
    Name: NotRequired[Literal["always", "on-failure"]]

@type_check_only
class _TopResult(TypedDict):
    Titles: list[str]
    Processes: list[list[str]]

class Container(Model):
    @property
    def name(self) -> str | None: ...
    @property
    def image(self) -> Image | None: ...
    @property
    def labels(self): ...
    @property
    def status(self) -> str: ...
    @property
    def health(self) -> str: ...
    @property
    def ports(self) -> dict[str, list[dict[str, str]] | None]: ...
    @overload
    def attach(
        self,
        *,
        stdout: bool = True,
        stderr: bool = True,
        stream: Literal[False] = False,
        logs: bool = False,
        demux: Literal[False] = False,
    ) -> bytes: ...
    @overload
    def attach(
        self,
        *,
        stdout: bool = True,
        stderr: bool = True,
        stream: Literal[False] = False,
        logs: bool = False,
        demux: Literal[True],
    ) -> tuple[bytes | None, bytes | None]: ...
    @overload
    def attach(
        self,
        *,
        stdout: bool = True,
        stderr: bool = True,
        stream: Literal[True],
        logs: bool = False,
        demux: Literal[False] = False,
    ) -> CancellableStream[bytes]: ...
    @overload
    def attach(
        self, *, stdout: bool = True, stderr: bool = True, stream: Literal[True], logs: bool = False, demux: Literal[True]
    ) -> CancellableStream[tuple[bytes | None, bytes | None]]: ...
    def attach_socket(self, **kwargs) -> SocketIO | _BufferedReaderStream | SSHSocket: ...
    def commit(self, repository: str | None = None, tag: str | None = None, **kwargs) -> Image: ...
    def diff(self) -> list[dict[str, int | str]]: ...
    def exec_run(
        self,
        cmd: str | list[str],
        stdout: bool = True,
        stderr: bool = True,
        stdin: bool = False,
        tty: bool = False,
        privileged: bool = False,
        user: str = "",
        detach: bool = False,
        stream: bool = False,
        socket: bool = False,
        environment: dict[str, str] | list[str] | None = None,
        workdir: str | None = None,
        demux: bool = False,
    ) -> ExecResult: ...
    def export(self, chunk_size: int | None = 2097152) -> str: ...
    def get_archive(
        self, path: str, chunk_size: int | None = 2097152, encode_stream: bool = False
    ) -> tuple[Iterator[bytes], dict[str, Any] | None]: ...
    def kill(self, signal: str | int | None = None) -> None: ...
    @overload
    def logs(
        self,
        *,
        stdout: bool = True,
        stderr: bool = True,
        stream: Literal[True],
        timestamps: bool = False,
        tail: Literal["all"] | int = "all",
        since: datetime.datetime | float | None = None,
        follow: bool | None = None,
        until: datetime.datetime | float | None = None,
    ) -> CancellableStream[bytes]: ...
    @overload
    def logs(
        self,
        *,
        stdout: bool = True,
        stderr: bool = True,
        stream: Literal[False] = False,
        timestamps: bool = False,
        tail: Literal["all"] | int = "all",
        since: datetime.datetime | float | None = None,
        follow: bool | None = None,
        until: datetime.datetime | float | None = None,
    ) -> bytes: ...
    def pause(self) -> None: ...
    def put_archive(self, path: str, data) -> bool: ...
    def remove(self, *, v: bool = False, link: bool = False, force: bool = False) -> None: ...
    def rename(self, name: str) -> None: ...
    def resize(self, height: int, width: int) -> None: ...
    def restart(self, *, timeout: float | None = 10) -> None: ...
    def start(self) -> None: ...
    def stats(self, **kwargs) -> Iterator[dict[str, Any]] | dict[str, Any]: ...
    def stop(self, *, timeout: float | None = None) -> None: ...
    def top(self, *, ps_args: str | None = None) -> _TopResult: ...
    def unpause(self) -> None: ...
    def update(
        self,
        *,
        blkio_weight: int | None = None,
        cpu_period: int | None = None,
        cpu_quota: int | None = None,
        cpu_shares: int | None = None,
        cpuset_cpus: str | None = None,
        cpuset_mems: str | None = None,
        mem_limit: float | str | None = None,
        mem_reservation: float | str | None = None,
        memswap_limit: int | str | None = None,
        kernel_memory: int | str | None = None,
        restart_policy: _RestartPolicy | None = None,
    ): ...
    def wait(
        self, *, timeout: float | None = None, condition: Literal["not-running", "next-exit", "removed"] | None = None
    ) -> WaitContainerResponse: ...

class ContainerCollection(Collection[Container]):
    model: type[Container]
    @overload
    def run(
        self,
        image: str | Image,
        command: str | list[str] | None = None,
        stdout: bool = True,
        stderr: bool = False,
        remove: bool = False,
        *,
        auto_remove: bool = False,
        blkio_weight_device: list[ContainerWeightDevice] | None = None,
        blkio_weight: int | None = None,
        cap_add: list[str] | None = None,
        cap_drop: list[str] | None = None,
        cgroup_parent: str | None = None,
        cgroupns: Literal["private", "host"] | None = None,
        cpu_count: int | None = None,
        cpu_percent: int | None = None,
        cpu_period: int | None = None,
        cpu_quota: int | None = None,
        cpu_rt_period: int | None = None,
        cpu_rt_runtime: int | None = None,
        cpu_shares: int | None = None,
        cpuset_cpus: str | None = None,
        cpuset_mems: str | None = None,
        detach: Literal[False] = False,
        device_cgroup_rules: list[str] | None = None,
        device_read_bps: list[Mapping[str, str | int]] | None = None,
        device_read_iops: list[Mapping[str, str | int]] | None = None,
        device_write_bps: list[Mapping[str, str | int]] | None = None,
        device_write_iops: list[Mapping[str, str | int]] | None = None,
        devices: list[str] | None = None,
        device_requests: list[DeviceRequest] | None = None,
        dns: list[str] | None = None,
        dns_opt: list[str] | None = None,
        dns_search: list[str] | None = None,
        domainname: str | list[str] | None = None,
        entrypoint: str | list[str] | None = None,
        environment: dict[str, str] | list[str] | None = None,
        extra_hosts: dict[str, str] | None = None,
        group_add: Iterable[str | int] | None = None,
        healthcheck: dict[str, Any] | None = None,
        hostname: str | None = None,
        init: bool | None = None,
        init_path: str | None = None,
        ipc_mode: str | None = None,
        isolation: str | None = None,
        kernel_memory: str | int | None = None,
        labels: dict[str, str] | list[str] | None = None,
        links: dict[str, str] | dict[str, None] | dict[str, str | None] | Iterable[tuple[str, str | None]] | None = None,
        log_config: LogConfig | None = None,
        lxc_conf: dict[str, str] | None = None,
        mac_address: str | None = None,
        mem_limit: str | int | None = None,
        mem_reservation: str | int | None = None,
        mem_swappiness: int | None = None,
        memswap_limit: str | int | None = None,
        mounts: list[Mount] | None = None,
        name: str | None = None,
        nano_cpus: int | None = None,
        network: str | None = None,
        network_disabled: bool = False,
        network_mode: str | None = None,
        networking_config: dict[str, EndpointConfig] | None = None,
        oom_kill_disable: bool = False,
        oom_score_adj: int | None = None,
        pid_mode: str | None = None,
        pids_limit: int | None = None,
        platform: str | None = None,
        ports: Mapping[str, int | list[int] | tuple[str, int] | None] | None = None,
        privileged: bool = False,
        publish_all_ports: bool = False,
        read_only: bool | None = None,
        restart_policy: _RestartPolicy | None = None,
        runtime: str | None = None,
        security_opt: list[str] | None = None,
        shm_size: str | int | None = None,
        stdin_open: bool = False,
        stop_signal: str | None = None,
        storage_opt: dict[str, str] | None = None,
        stream: bool = False,
        sysctls: dict[str, str] | None = None,
        tmpfs: dict[str, str] | None = None,
        tty: bool = False,
        ulimits: list[Ulimit] | None = None,
        use_config_proxy: bool | None = None,
        user: str | int | None = None,
        userns_mode: str | None = None,
        uts_mode: str | None = None,
        version: str | None = None,
        volume_driver: str | None = None,
        volumes: dict[str, dict[str, str]] | list[str] | None = None,
        volumes_from: list[str] | None = None,
        working_dir: str | None = None,
    ) -> bytes: ...  # TODO: This should return a stream, if `stream` is True
    @overload
    def run(
        self,
        image: str | Image,
        command: str | list[str] | None = None,
        stdout: bool = True,
        stderr: bool = False,
        remove: bool = False,
        *,
        auto_remove: bool = False,
        blkio_weight_device: list[ContainerWeightDevice] | None = None,
        blkio_weight: int | None = None,
        cap_add: list[str] | None = None,
        cap_drop: list[str] | None = None,
        cgroup_parent: str | None = None,
        cgroupns: Literal["private", "host"] | None = None,
        cpu_count: int | None = None,
        cpu_percent: int | None = None,
        cpu_period: int | None = None,
        cpu_quota: int | None = None,
        cpu_rt_period: int | None = None,
        cpu_rt_runtime: int | None = None,
        cpu_shares: int | None = None,
        cpuset_cpus: str | None = None,
        cpuset_mems: str | None = None,
        detach: Literal[True],
        device_cgroup_rules: list[str] | None = None,
        device_read_bps: list[Mapping[str, str | int]] | None = None,
        device_read_iops: list[Mapping[str, str | int]] | None = None,
        device_write_bps: list[Mapping[str, str | int]] | None = None,
        device_write_iops: list[Mapping[str, str | int]] | None = None,
        devices: list[str] | None = None,
        device_requests: list[DeviceRequest] | None = None,
        dns: list[str] | None = None,
        dns_opt: list[str] | None = None,
        dns_search: list[str] | None = None,
        domainname: str | list[str] | None = None,
        entrypoint: str | list[str] | None = None,
        environment: dict[str, str] | list[str] | None = None,
        extra_hosts: dict[str, str] | None = None,
        group_add: Iterable[str | int] | None = None,
        healthcheck: dict[str, Any] | None = None,
        hostname: str | None = None,
        init: bool | None = None,
        init_path: str | None = None,
        ipc_mode: str | None = None,
        isolation: str | None = None,
        kernel_memory: str | int | None = None,
        labels: dict[str, str] | list[str] | None = None,
        links: dict[str, str] | dict[str, None] | dict[str, str | None] | Iterable[tuple[str, str | None]] | None = None,
        log_config: LogConfig | None = None,
        lxc_conf: dict[str, str] | None = None,
        mac_address: str | None = None,
        mem_limit: str | int | None = None,
        mem_reservation: str | int | None = None,
        mem_swappiness: int | None = None,
        memswap_limit: str | int | None = None,
        mounts: list[Mount] | None = None,
        name: str | None = None,
        nano_cpus: int | None = None,
        network: str | None = None,
        network_disabled: bool = False,
        network_mode: str | None = None,
        networking_config: dict[str, EndpointConfig] | None = None,
        oom_kill_disable: bool = False,
        oom_score_adj: int | None = None,
        pid_mode: str | None = None,
        pids_limit: int | None = None,
        platform: str | None = None,
        ports: Mapping[str, int | list[int] | tuple[str, int] | None] | None = None,
        privileged: bool = False,
        publish_all_ports: bool = False,
        read_only: bool | None = None,
        restart_policy: _RestartPolicy | None = None,
        runtime: str | None = None,
        security_opt: list[str] | None = None,
        shm_size: str | int | None = None,
        stdin_open: bool = False,
        stop_signal: str | None = None,
        storage_opt: dict[str, str] | None = None,
        stream: bool = False,
        sysctls: dict[str, str] | None = None,
        tmpfs: dict[str, str] | None = None,
        tty: bool = False,
        ulimits: list[Ulimit] | None = None,
        use_config_proxy: bool | None = None,
        user: str | int | None = None,
        userns_mode: str | None = None,
        uts_mode: str | None = None,
        version: str | None = None,
        volume_driver: str | None = None,
        volumes: dict[str, dict[str, str]] | list[str] | None = None,
        volumes_from: list[str] | None = None,
        working_dir: str | None = None,
    ) -> Container: ...
    @override
    def create(  # type: ignore[override]
        self,
        image: str | Image,
        command: str | list[str] | None = None,
        *,
        auto_remove: bool = False,
        blkio_weight_device: list[ContainerWeightDevice] | None = None,
        blkio_weight: int | None = None,
        cap_add: list[str] | None = None,
        cap_drop: list[str] | None = None,
        cgroup_parent: str | None = None,
        cgroupns: Literal["private", "host"] | None = None,
        cpu_count: int | None = None,
        cpu_percent: int | None = None,
        cpu_period: int | None = None,
        cpu_quota: int | None = None,
        cpu_rt_period: int | None = None,
        cpu_rt_runtime: int | None = None,
        cpu_shares: int | None = None,
        cpuset_cpus: str | None = None,
        cpuset_mems: str | None = None,
        detach: bool = False,
        device_cgroup_rules: list[str] | None = None,
        device_read_bps: list[Mapping[str, str | int]] | None = None,
        device_read_iops: list[Mapping[str, str | int]] | None = None,
        device_write_bps: list[Mapping[str, str | int]] | None = None,
        device_write_iops: list[Mapping[str, str | int]] | None = None,
        devices: list[str] | None = None,
        device_requests: list[DeviceRequest] | None = None,
        dns: list[str] | None = None,
        dns_opt: list[str] | None = None,
        dns_search: list[str] | None = None,
        domainname: str | list[str] | None = None,
        entrypoint: str | list[str] | None = None,
        environment: dict[str, str] | list[str] | None = None,
        extra_hosts: dict[str, str] | None = None,
        group_add: Iterable[str | int] | None = None,
        healthcheck: dict[str, Any] | None = None,
        hostname: str | None = None,
        init: bool | None = None,
        init_path: str | None = None,
        ipc_mode: str | None = None,
        isolation: str | None = None,
        kernel_memory: str | int | None = None,
        labels: dict[str, str] | list[str] | None = None,
        links: dict[str, str] | dict[str, None] | dict[str, str | None] | Iterable[tuple[str, str | None]] | None = None,
        log_config: LogConfig | None = None,
        lxc_conf: dict[str, str] | None = None,
        mac_address: str | None = None,
        mem_limit: str | int | None = None,
        mem_reservation: str | int | None = None,
        mem_swappiness: int | None = None,
        memswap_limit: str | int | None = None,
        mounts: list[Mount] | None = None,
        name: str | None = None,
        nano_cpus: int | None = None,
        network: str | None = None,
        network_disabled: bool = False,
        network_mode: str | None = None,
        networking_config: dict[str, EndpointConfig] | None = None,
        oom_kill_disable: bool = False,
        oom_score_adj: int | None = None,
        pid_mode: str | None = None,
        pids_limit: int | None = None,
        platform: str | None = None,
        ports: Mapping[str, int | list[int] | tuple[str, int] | None] | None = None,
        privileged: bool = False,
        publish_all_ports: bool = False,
        read_only: bool | None = None,
        restart_policy: _RestartPolicy | None = None,
        runtime: str | None = None,
        security_opt: list[str] | None = None,
        shm_size: str | int | None = None,
        stdin_open: bool = False,
        stop_signal: str | None = None,
        storage_opt: dict[str, str] | None = None,
        stream: bool = False,
        sysctls: dict[str, str] | None = None,
        tmpfs: dict[str, str] | None = None,
        tty: bool = False,
        ulimits: list[Ulimit] | None = None,
        use_config_proxy: bool | None = None,
        user: str | int | None = None,
        userns_mode: str | None = None,
        uts_mode: str | None = None,
        version: str | None = None,
        volume_driver: str | None = None,
        volumes: dict[str, dict[str, str]] | list[str] | None = None,
        volumes_from: list[str] | None = None,
        working_dir: str | None = None,
    ) -> Container: ...
    @override
    def get(self, container_id: str) -> Container: ...
    def list(
        self,
        all: bool = False,
        before: str | None = None,
        filters: dict[str, str | list[str] | bool] | None = None,
        limit: int = -1,
        since: str | None = None,
        sparse: bool = False,
        ignore_removed: bool = False,
    ) -> list[Container]: ...
    def prune(self, filters: dict[str, Any] | None = None) -> dict[str, Any]: ...

RUN_CREATE_KWARGS: list[str]
RUN_HOST_CONFIG_KWARGS: list[str]

class ExecResult(NamedTuple):
    exit_code: int | None
    output: bytes | Iterator[bytes]
