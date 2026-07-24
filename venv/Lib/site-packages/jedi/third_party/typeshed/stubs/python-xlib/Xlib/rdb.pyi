from _typeshed import SupportsDunderGT, SupportsDunderLT, SupportsRead
from collections.abc import Iterable, Mapping, Sequence
from re import Pattern
from typing import Any, Final, Protocol, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from Xlib.display import Display
from Xlib.support.lock import _DummyLock

_T = TypeVar("_T")
_T_contra = TypeVar("_T_contra", contravariant=True)

_DB: TypeAlias = dict[str, tuple[_DB, ...]]
# A recursive type can be a bit annoying due to dict invariance,
# so this is a slightly less precise version of the _DB alias for parameter annotations
_DB_Param: TypeAlias = dict[str, Any]

@type_check_only
class _SupportsComparisons(SupportsDunderLT[_T_contra], SupportsDunderGT[_T_contra], Protocol[_T_contra]): ...

comment_re: Final[Pattern[str]]
resource_spec_re: Final[Pattern[str]]
value_escape_re: Final[Pattern[str]]
resource_parts_re: Final[Pattern[str]]
NAME_MATCH: Final = 0
CLASS_MATCH: Final = 2
WILD_MATCH: Final = 4
MATCH_SKIP: Final = 6

class OptionError(Exception): ...

class ResourceDB:
    db: _DB
    lock: _DummyLock
    def __init__(
        self,
        file: bytes | SupportsRead[str] | None = None,
        string: str | None = None,
        resources: Iterable[tuple[str, object]] | None = None,
    ) -> None: ...
    def insert_file(self, file: bytes | SupportsRead[str]) -> None: ...
    def insert_string(self, data: str) -> None: ...
    def insert_resources(self, resources: Iterable[tuple[str, object]]) -> None: ...
    def insert(self, resource: str, value: object) -> None: ...
    def __getitem__(self, keys_tuple: tuple[str, str]) -> Any: ...
    @overload
    def get(self, res: str, cls: str, default: None = None) -> Any: ...
    @overload
    def get(self, res: str, cls: str, default: _T) -> _T: ...
    def update(self, db: ResourceDB) -> None: ...
    def output(self) -> str: ...
    def getopt(self, name: str, argv: Sequence[str], opts: Mapping[str, Option]) -> Sequence[str]: ...

def bin_insert(list: list[_SupportsComparisons[_T]], element: _SupportsComparisons[_T]) -> None: ...
def update_db(dest: _DB_Param, src: _DB_Param) -> None: ...
def copy_group(group: tuple[_DB_Param, ...]) -> tuple[_DB, ...]: ...
def copy_db(db: _DB_Param) -> _DB: ...
def output_db(prefix: str, db: _DB_Param) -> str: ...
def output_escape(value: object) -> str: ...

class Option:
    def parse(self, name: str, db: ResourceDB, args: Sequence[_T]) -> Sequence[_T]: ...

class NoArg(Option):
    specifier: str
    value: object
    def __init__(self, specifier: str, value: object) -> None: ...

class IsArg(Option):
    specifier: str
    def __init__(self, specifier: str) -> None: ...

class SepArg(Option):
    specifier: str
    def __init__(self, specifier: str) -> None: ...

class ResArgClass(Option):
    def parse(self, name: str, db: ResourceDB, args: Sequence[str]) -> Sequence[str]: ...  # type: ignore[override]

ResArg: ResArgClass

class SkipArgClass(Option): ...

SkipArg: SkipArgClass

class SkipLineClass(Option): ...

SkipLine: SkipLineClass

class SkipNArgs(Option):
    count: int
    def __init__(self, count: int) -> None: ...

def get_display_opts(
    options: Mapping[str, Option], argv: Sequence[str] = ...
) -> tuple[Display, str, ResourceDB, Sequence[str]]: ...

stdopts: Final[dict[str, SepArg | NoArg | ResArgClass]]
