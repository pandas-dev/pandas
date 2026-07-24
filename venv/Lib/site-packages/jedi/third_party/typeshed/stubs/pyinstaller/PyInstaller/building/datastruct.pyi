# https://pyinstaller.org/en/stable/advanced-topics.html#the-toc-and-tree-classes
from collections.abc import Iterable, Sequence
from typing import ClassVar, Literal, SupportsIndex
from typing_extensions import LiteralString, Self, TypeAlias

_TypeCode: TypeAlias = Literal["DEPENDENCY", "SYMLINK", "DATA", "BINARY", "EXECUTABLE", "EXTENSION", "OPTION"]
_TOCTuple: TypeAlias = tuple[str, str | None, _TypeCode | None]

class TOC(list[_TOCTuple]):
    filenames: set[str]
    def __init__(self, initlist: Iterable[_TOCTuple] | None = None) -> None: ...
    def append(self, entry: _TOCTuple) -> None: ...
    def insert(self, pos: SupportsIndex, entry: _TOCTuple) -> None: ...
    def __add__(self, other: Iterable[_TOCTuple]) -> TOC: ...  # type: ignore[override]
    def __radd__(self, other: Iterable[_TOCTuple]) -> TOC: ...
    def __iadd__(self, other: Iterable[_TOCTuple]) -> Self: ...  # type: ignore[override]
    def extend(self, other: Iterable[_TOCTuple]) -> None: ...
    def __sub__(self, other: Iterable[_TOCTuple]) -> TOC: ...
    def __rsub__(self, other: Iterable[_TOCTuple]) -> TOC: ...
    # slicing a TOC is not supported, but has a special case for slice(None, None, None)
    def __setitem__(self, key: int | slice, value: Iterable[_TOCTuple]) -> None: ...  # type: ignore[override]

class Target:
    invcnum: ClassVar[int]
    tocfilename: LiteralString
    tocbasename: LiteralString
    dependencies: list[_TOCTuple]
    def __init__(self) -> None: ...
    def __postinit__(self) -> None: ...

class Tree(Target, list[_TOCTuple]):
    root: str | None
    prefix: str | None
    excludes: Sequence[str]
    typecode: _TypeCode
    def __init__(
        self,
        root: str | None = None,
        prefix: str | None = None,
        excludes: Sequence[str] | None = None,
        typecode: _TypeCode = "DATA",
    ) -> None: ...
    def assemble(self) -> None: ...

def normalize_toc(toc: Iterable[_TOCTuple]) -> list[_TOCTuple]: ...
def normalize_pyz_toc(toc: Iterable[_TOCTuple]) -> list[_TOCTuple]: ...
def toc_process_symbolic_links(toc: Iterable[_TOCTuple]) -> list[_TOCTuple]: ...
