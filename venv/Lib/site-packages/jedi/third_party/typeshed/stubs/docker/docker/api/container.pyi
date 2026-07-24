import datetime
from _typeshed import Incomplete
from typing import Any, Literal, TypedDict, overload, type_check_only
from typing_extensions import TypeAlias

from docker._types import WaitContainerResponse
from docker.types.daemon import CancellableStream

from ..types import ContainerConfig, EndpointConfig, HostConfig, NetworkingConfig

@type_check_only
class _HasId(TypedDict):
    Id: str

@type_check_only
class _HasID(TypedDict):
    ID: str

@type_check_only
class _TopResult(TypedDict):
    Titles: list[str]
    Processes: list[list[str]]

_Container: TypeAlias = _HasId | _HasID | str

class ContainerApiMixin:
    @overload
    def attach(
        self,
        container: _Container,
        stdout: bool = True,
        stderr: bool = True,
        stream: Literal[False] = False,
        logs: bool = False,
        demux: Literal[False] = False,
    ) -> bytes: ...
    @overload
    def attach(
        self,
        container: _Container,
        stdout: bool = True,
        stderr: bool = True,
        stream: Literal[False] = False,
        logs: bool = False,
        *,
        demux: Literal[True],
    ) -> tuple[bytes | None, bytes | None]: ...
    @overload
    def attach(
        self,
        container: _Container,
        stdout: bool = True,
        stderr: bool = True,
        *,
        stream: Literal[True],
        logs: bool = False,
        demux: Literal[False] = False,
    ) -> CancellableStream[bytes]: ...
    @overload
    def attach(
        self,
        container: _Container,
        stdout: bool = True,
        stderr: bool = True,
        *,
        stream: Literal[True],
        logs: bool = False,
        demux: Literal[True],
    ) -> CancellableStream[tuple[bytes | None, bytes | None]]: ...
    def attach_socket(self, container: _Container, params=None, ws: bool = False): ...
    def commit(
        self,
        container: _Container,
        repository: str | None = None,
        tag: str | None = None,
        message=None,
        author=None,
        pause: bool = True,
        changes=None,
        conf=None,
    ): ...
    def containers(
        self,
        quiet: bool = False,
        all: bool = False,
        trunc: bool = False,
        latest: bool = False,
        since: str | None = None,
        before: str | None = None,
        limit: int = -1,
        size: bool = False,
        filters=None,
    ): ...
    def create_container(
        self,
        image,
        command: str | list[str] | None = None,
        hostname: str | None = None,
        user: str | int | None = None,
        detach: bool = False,
        stdin_open: bool = False,
        tty: bool = False,
        # list is invariant, enumerating all possible union combination would be too complex for:
        # list[str | int | tuple[int | str, str] | tuple[int | str, ...]]
        ports: dict[str, dict[Incomplete, Incomplete]] | list[Any] | None = None,
        environment: dict[str, str] | list[str] | None = None,
        volumes: str | list[str] | None = None,
        network_disabled: bool = False,
        name: str | None = None,
        entrypoint: str | list[str] | None = None,
        working_dir: str | None = None,
        domainname: str | None = None,
        host_config=None,
        mac_address: str | None = None,
        labels: dict[str, str] | list[str] | None = None,
        stop_signal: str | None = None,
        networking_config=None,
        healthcheck=None,
        stop_timeout: int | None = None,
        runtime: str | None = None,
        use_config_proxy: bool = True,
        platform: str | None = None,
    ): ...
    def create_container_config(self, *args, **kwargs) -> ContainerConfig: ...
    def create_container_from_config(self, config, name=None, platform=None): ...
    def create_host_config(self, *args, **kwargs) -> HostConfig: ...
    def create_networking_config(self, *args, **kwargs) -> NetworkingConfig: ...
    def create_endpoint_config(self, *args, **kwargs) -> EndpointConfig: ...
    def diff(self, container: _Container) -> list[dict[Incomplete, Incomplete]]: ...
    def export(self, container: _Container, chunk_size: int | None = 2097152): ...
    def get_archive(
        self, container: _Container, path, chunk_size: int | None = 2097152, encode_stream: bool = False
    ) -> tuple[Incomplete, Incomplete]: ...
    def inspect_container(self, container: _Container): ...
    def kill(self, container: _Container, signal: str | int | None = None) -> None: ...
    @overload
    def logs(
        self,
        container: _Container,
        stdout: bool = True,
        stderr: bool = True,
        *,
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
        container: _Container,
        stdout: bool,
        stderr: bool,
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
        container: _Container,
        stdout: bool = True,
        stderr: bool = True,
        stream: Literal[False] = False,
        timestamps: bool = False,
        tail: Literal["all"] | int = "all",
        since: datetime.datetime | float | None = None,
        follow: bool | None = None,
        until: datetime.datetime | float | None = None,
    ) -> bytes: ...
    def pause(self, container: _Container) -> None: ...
    def port(self, container: _Container, private_port: int): ...
    def put_archive(self, container: _Container, path: str, data) -> bool: ...
    def prune_containers(self, filters=None): ...
    def remove_container(self, container: _Container, v: bool = False, link: bool = False, force: bool = False) -> None: ...
    def rename(self, container: _Container, name: str) -> None: ...
    def resize(self, container: _Container, height: int, width: int) -> None: ...
    def restart(self, container: _Container, timeout: int = 10) -> None: ...
    def start(self, container: _Container) -> None: ...
    def stats(self, container: _Container, decode: bool | None = None, stream: bool = True, one_shot: bool | None = None): ...
    def stop(self, container: _Container, timeout: int | None = None) -> None: ...
    def top(self, container: _Container, ps_args: str | None = None) -> _TopResult: ...
    def unpause(self, container: _Container) -> None: ...
    def update_container(
        self,
        container: _Container,
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
        restart_policy=None,
    ): ...
    def wait(
        self,
        container: _Container,
        timeout: int | None = None,
        condition: Literal["not-running", "next-exit", "removed"] | None = None,
    ) -> WaitContainerResponse: ...
