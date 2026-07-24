from _typeshed import SupportsRead
from collections.abc import Iterator
from io import StringIO
from typing import IO, Any, Literal, TypedDict, overload, type_check_only
from typing_extensions import TypeAlias

from docker._types import JSON

from .resource import Collection, Model

_ImageList: TypeAlias = list[Image]  # To resolve conflicts with a method called "list"

@type_check_only
class _ContainerLimits(TypedDict, total=False):
    memory: int
    memswap: int
    cpushares: int
    cpusetcpus: str

class Image(Model):
    @property
    def labels(self) -> dict[str, Any]: ...
    @property
    def short_id(self) -> str: ...
    @property
    def tags(self) -> list[str]: ...
    def history(self) -> list[Any]: ...
    def remove(self, force: bool = False, noprune: bool = False) -> dict[str, Any]: ...
    def save(self, chunk_size: int = 2097152, named: str | bool = False) -> Iterator[Any]: ...
    def tag(self, repository: str, tag: str | None = None, **kwargs) -> bool: ...

class RegistryData(Model):
    image_name: str
    def __init__(self, image_name: str, *args, **kwargs) -> None: ...
    @property
    def id(self) -> str: ...
    @property
    def short_id(self) -> str: ...
    def pull(self, platform: str | None = None) -> Image: ...
    def has_platform(self, platform): ...
    def reload(self) -> None: ...

class ImageCollection(Collection[Image]):
    model: type[Image]
    def build(
        self,
        *,
        path: str | None = None,
        fileobj: StringIO | IO[bytes] | None = None,
        tag: str | None = None,
        quiet: bool = False,
        nocache: bool = False,
        rm: bool = False,
        timeout: int | None = None,
        custom_context: bool = False,
        encoding: str | None = None,
        pull: bool = False,
        forcerm: bool = False,
        dockerfile: str | None = None,
        buildargs: dict[str, Any] | None = None,
        container_limits: _ContainerLimits | None = None,
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
    ) -> tuple[Image, Iterator[JSON]]: ...
    def get(self, name: str) -> Image: ...
    def get_registry_data(self, name, auth_config: dict[str, Any] | None = None) -> RegistryData: ...
    def list(self, name: str | None = None, all: bool = False, filters: dict[str, Any] | None = None) -> _ImageList: ...
    def load(self, data: bytes | SupportsRead[bytes]) -> _ImageList: ...
    @overload
    def pull(
        self,
        repository: str,
        tag: str | None = None,
        all_tags: Literal[False] = False,
        *,
        platform: str | None = None,
        auth_config: dict[str, Any] | None = None,
    ) -> Image: ...
    @overload
    def pull(
        self,
        repository: str,
        tag: str | None = None,
        *,
        all_tags: Literal[True],
        auth_config: dict[str, Any] | None = None,
        platform: str | None = None,
    ) -> _ImageList: ...
    @overload
    def pull(
        self,
        repository: str,
        tag: str | None,
        all_tags: Literal[True],
        *,
        auth_config: dict[str, Any] | None = None,
        platform: str | None = None,
    ) -> _ImageList: ...
    def push(self, repository: str, tag: str | None = None, **kwargs): ...
    def remove(self, *args, **kwargs) -> None: ...
    def search(self, *args, **kwargs): ...
    def prune(self, filters: dict[str, Any] | None = None): ...
    def prune_builds(self, *args, **kwargs): ...

def normalize_platform(platform, engine_info): ...
