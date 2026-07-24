from collections.abc import Generator
from io import StringIO
from logging import Logger
from typing import IO, Any, TypedDict, type_check_only

log: Logger

@type_check_only
class _ContainerLimits(TypedDict, total=False):
    memory: int
    memswap: int
    cpushares: int
    cpusetcpus: str

@type_check_only
class _Filers(TypedDict, total=False):
    dangling: bool
    until: str

class BuildApiMixin:
    def build(
        self,
        path: str | None = None,
        tag: str | None = None,
        quiet: bool = False,
        fileobj: StringIO | IO[bytes] | None = None,
        nocache: bool = False,
        rm: bool = False,
        timeout: int | None = None,
        custom_context: bool = False,
        encoding: str | None = None,
        pull: bool = False,
        forcerm: bool = False,
        dockerfile: str | None = None,
        container_limits: _ContainerLimits | None = None,
        decode: bool = False,
        buildargs: dict[str, Any] | None = None,
        gzip: bool = False,
        shmsize: int | None = None,
        labels: dict[str, Any] | None = None,
        # need to use list, because the type must be json serializable
        cache_from: list[str] | None = None,
        target: str | None = None,
        network_mode: str | None = None,
        squash: bool | None = None,
        extra_hosts: list[str] | dict[str, str] | None = None,
        platform: str | None = None,
        isolation: str | None = None,
        use_config_proxy: bool = True,
    ) -> Generator[Any]: ...
    def prune_builds(
        self, filters: _Filers | None = None, keep_storage: int | None = None, all: bool | None = None
    ) -> dict[str, Any]: ...

def process_dockerfile(dockerfile: str | None, path: str) -> tuple[str | None, str | None]: ...
