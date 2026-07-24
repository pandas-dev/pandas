from _typeshed import Incomplete
from collections.abc import Generator, Iterable
from os import PathLike
from typing import TypeVar
from typing_extensions import TypeAlias

_DirNode: TypeAlias = Incomplete  # scantree.DirNode
_RecursionPath: TypeAlias = Incomplete  # scantree.RecursionPath
_RP = TypeVar("_RP", bound=_RecursionPath)

__all__ = [
    "__version__",
    "algorithms_guaranteed",
    "algorithms_available",
    "dirhash",
    "dirhash_impl",
    "included_paths",
    "Filter",
    "get_match_patterns",
    "Protocol",
]

__version__: str
algorithms_guaranteed: set[str]
algorithms_available: set[str]

def dirhash(
    directory: str | PathLike[str],
    algorithm: str,
    match: Iterable[str] = ("*",),
    ignore: Iterable[str] | None = None,
    linked_dirs: bool = True,
    linked_files: bool = True,
    empty_dirs: bool = False,
    entry_properties: Iterable[str] = ("name", "data"),
    allow_cyclic_links: bool = False,
    chunk_size: int = 1048576,
    jobs: int = 1,
) -> str: ...
def dirhash_impl(
    directory: str | PathLike[str],
    algorithm: str,
    filter_: Filter | None = None,
    protocol: Protocol | None = None,
    chunk_size: int = 1048576,
    jobs: int = 1,
) -> str: ...
def included_paths(
    directory: str | PathLike[str],
    match: Iterable[str] = ("*",),
    ignore: Iterable[str] | None = None,
    linked_dirs: bool = True,
    linked_files: bool = True,
    empty_dirs: bool = False,
    allow_cyclic_links: bool = False,
) -> list[str]: ...

class Filter:
    linked_dirs: bool
    linked_files: bool
    empty_dirs: bool

    def __init__(
        self,
        match_patterns: Iterable[str] | None = None,
        linked_dirs: bool = True,
        linked_files: bool = True,
        empty_dirs: bool = False,
    ) -> None: ...
    @property
    def match_patterns(self) -> tuple[str, ...]: ...
    def include(self, recursion_path: _RecursionPath) -> bool: ...
    def match_file(self, filepath: str | PathLike[str]) -> bool: ...
    def __call__(self, paths: Iterable[_RP]) -> Generator[_RP]: ...

def get_match_patterns(
    match: Iterable[str] | None = None,
    ignore: Iterable[str] | None = None,
    ignore_extensions: Iterable[str] | None = None,
    ignore_hidden: bool = False,
) -> list[str]: ...

class Protocol:
    class EntryProperties:
        NAME: str
        DATA: str
        IS_LINK: str
        options: set[str]

    entry_properties: Iterable[str]
    allow_cyclic_links: bool
    def __init__(self, entry_properties: Iterable[str] = ("name", "data"), allow_cyclic_links: bool = False) -> None: ...
    def get_descriptor(self, dir_node: _DirNode) -> str: ...
