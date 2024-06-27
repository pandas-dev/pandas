from _typeshed import Incomplete
from collections.abc import Callable, Iterable
from types import TracebackType
from typing import Any, ClassVar, Final, Literal

from setuptools import Command

def safe_name(name: str) -> str: ...
def safe_version(version: str) -> str: ...

setuptools_major_version: Final[int]

PY_LIMITED_API_PATTERN: Final[str]

def python_tag() -> str: ...
def get_platform(archive_root: str | None) -> str: ...
def get_flag(var: str, fallback: bool, expected: bool = True, warn: bool = True) -> bool: ...
def get_abi_tag() -> str | None: ...
def safer_name(name: str) -> str: ...
def safer_version(version: str) -> str: ...
def remove_readonly(
    func: Callable[..., object], path: str, excinfo: tuple[type[Exception], Exception, TracebackType]
) -> None: ...
def remove_readonly_exc(func: Callable[..., object], path: str, exc: Exception) -> None: ...

class bdist_wheel(Command):
    description: ClassVar[str]
    supported_compressions: ClassVar[dict[str, int]]
    user_options: ClassVar[list[tuple[Any, ...]]]
    boolean_options: ClassVar[list[str]]

    bdist_dir: str | None
    data_dir: Incomplete | None
    plat_name: str | None
    plat_tag: Incomplete | None
    format: str
    keep_temp: bool
    dist_dir: str | None
    egginfo_dir: Incomplete | None
    root_is_pure: bool | None
    skip_build: Incomplete | None
    relative: bool
    owner: Incomplete | None
    group: Incomplete | None
    universal: bool
    compression: str | int
    python_tag: str
    build_number: str | None
    py_limited_api: str | Literal[False]
    plat_name_supplied: bool

    def initialize_options(self) -> None: ...
    def finalize_options(self) -> None: ...
    @property
    def wheel_dist_name(self) -> str: ...
    def get_tag(self) -> tuple[str, str, str]: ...
    def run(self) -> None: ...
    def write_wheelfile(self, wheelfile_base: str, generator: str = ...) -> None: ...
    @property
    def license_paths(self) -> Iterable[str]: ...
    def egg2dist(self, egginfo_path: str, distinfo_path: str) -> None: ...
