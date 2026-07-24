import logging
from collections.abc import Mapping, Sequence
from typing import IO, ClassVar, TypedDict, type_check_only

from .util import Context

logger: logging.Logger

class KeyValues(dict[str, str]):
    parser_attr: ClassVar[str | None]
    parser: DockerfileParser
    def __init__(self, key_values: Mapping[str, str], parser: DockerfileParser) -> None: ...
    def __delitem__(self, key: str) -> None: ...
    def __setitem__(self, key: str, value: str) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...  # type: ignore[override]

class Labels(KeyValues): ...
class Envs(KeyValues): ...
class Args(KeyValues): ...

@type_check_only
class _InstructionDict(TypedDict):
    instruction: str
    startline: int
    endline: int
    content: str
    value: str

class DockerfileParser:
    fileobj: IO[str]
    dockerfile_path: str
    cache_content: bool
    cached_content: str
    env_replace: bool
    parent_env: dict[str, str]
    build_args: dict[str, str]
    def __init__(
        self,
        path: str | None = None,
        cache_content: bool = False,
        env_replace: bool = True,
        parent_env: dict[str, str] | None = None,
        fileobj: IO[str] | None = None,
        build_args: dict[str, str] | None = None,
    ) -> None: ...
    lines: list[str]
    content: str
    @property
    def structure(self) -> list[_InstructionDict]: ...
    @property
    def json(self) -> str: ...
    parent_images: Sequence[str]
    @property
    def is_multistage(self) -> bool: ...
    baseimage: str
    cmd: str
    labels: Mapping[str, str]
    envs: Mapping[str, str]
    args: Mapping[str, str]
    def add_lines(
        self, *lines: str, all_stages: bool | None = ..., at_start: bool | None = ..., skip_scratch: bool | None = ...
    ) -> None: ...
    def add_lines_at(
        self, anchor: str | int | dict[str, int], *lines: str, replace: bool | None = ..., after: bool | None = ...
    ) -> None: ...
    @property
    def context_structure(self) -> list[Context]: ...

def image_from(from_value: str) -> tuple[str | None, str | None]: ...
