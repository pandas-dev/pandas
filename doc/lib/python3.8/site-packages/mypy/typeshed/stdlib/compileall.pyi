import sys
from _typeshed import StrPath
from typing import Any, Optional, Protocol

if sys.version_info >= (3, 7):
    from py_compile import PycInvalidationMode

class _SupportsSearch(Protocol):
    def search(self, string: str) -> Any: ...

if sys.version_info >= (3, 9):
    def compile_dir(
        dir: StrPath,
        maxlevels: Optional[int] = ...,
        ddir: Optional[StrPath] = ...,
        force: bool = ...,
        rx: Optional[_SupportsSearch] = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        workers: int = ...,
        invalidation_mode: Optional[PycInvalidationMode] = ...,
        *,
        stripdir: Optional[str] = ...,  # TODO: change to Optional[StrPath] once https://bugs.python.org/issue40447 is resolved
        prependdir: Optional[StrPath] = ...,
        limit_sl_dest: Optional[StrPath] = ...,
        hardlink_dupes: bool = ...,
    ) -> int: ...
    def compile_file(
        fullname: StrPath,
        ddir: Optional[StrPath] = ...,
        force: bool = ...,
        rx: Optional[_SupportsSearch] = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        invalidation_mode: Optional[PycInvalidationMode] = ...,
        *,
        stripdir: Optional[str] = ...,  # TODO: change to Optional[StrPath] once https://bugs.python.org/issue40447 is resolved
        prependdir: Optional[StrPath] = ...,
        limit_sl_dest: Optional[StrPath] = ...,
        hardlink_dupes: bool = ...,
    ) -> int: ...

elif sys.version_info >= (3, 7):
    def compile_dir(
        dir: StrPath,
        maxlevels: int = ...,
        ddir: Optional[StrPath] = ...,
        force: bool = ...,
        rx: Optional[_SupportsSearch] = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        workers: int = ...,
        invalidation_mode: Optional[PycInvalidationMode] = ...,
    ) -> int: ...
    def compile_file(
        fullname: StrPath,
        ddir: Optional[StrPath] = ...,
        force: bool = ...,
        rx: Optional[_SupportsSearch] = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        invalidation_mode: Optional[PycInvalidationMode] = ...,
    ) -> int: ...

else:
    def compile_dir(
        dir: StrPath,
        maxlevels: int = ...,
        ddir: Optional[StrPath] = ...,
        force: bool = ...,
        rx: Optional[_SupportsSearch] = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        workers: int = ...,
    ) -> int: ...
    def compile_file(
        fullname: StrPath,
        ddir: Optional[StrPath] = ...,
        force: bool = ...,
        rx: Optional[_SupportsSearch] = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
    ) -> int: ...

if sys.version_info >= (3, 7):
    def compile_path(
        skip_curdir: bool = ...,
        maxlevels: int = ...,
        force: bool = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
        invalidation_mode: Optional[PycInvalidationMode] = ...,
    ) -> int: ...

else:
    def compile_path(
        skip_curdir: bool = ...,
        maxlevels: int = ...,
        force: bool = ...,
        quiet: int = ...,
        legacy: bool = ...,
        optimize: int = ...,
    ) -> int: ...
