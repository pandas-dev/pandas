from _typeshed import Incomplete, SupportsRead
from collections.abc import Mapping, Sequence
from typing import IO, Any
from typing_extensions import TypeAlias

from ._yaml import CEmitter, CParser
from .constructor import BaseConstructor, FullConstructor, SafeConstructor, UnsafeConstructor
from .representer import BaseRepresenter, SafeRepresenter
from .resolver import BaseResolver, Resolver

__all__ = ["CBaseLoader", "CSafeLoader", "CFullLoader", "CUnsafeLoader", "CLoader", "CBaseDumper", "CSafeDumper", "CDumper"]

_Readable: TypeAlias = SupportsRead[str | bytes]
_CLoader: TypeAlias = CLoader | CBaseLoader | CFullLoader | CSafeLoader | CUnsafeLoader  # noqa: Y047  # Used in other modules

class CBaseLoader(CParser, BaseConstructor, BaseResolver):
    def __init__(self, stream: str | bytes | _Readable) -> None: ...

class CLoader(CParser, SafeConstructor, Resolver):
    def __init__(self, stream: str | bytes | _Readable) -> None: ...

class CSafeLoader(CParser, SafeConstructor, Resolver):
    def __init__(self, stream: str | bytes | _Readable) -> None: ...

class CFullLoader(CParser, FullConstructor, Resolver):
    def __init__(self, stream: str | bytes | _Readable) -> None: ...

class CUnsafeLoader(CParser, UnsafeConstructor, Resolver):
    def __init__(self, stream: str | bytes | _Readable) -> None: ...

class CBaseDumper(CEmitter, BaseRepresenter, BaseResolver):
    def __init__(
        self,
        stream: IO[Any],
        default_style: str | None = None,
        default_flow_style: bool | None = False,
        canonical: Incomplete | None = None,
        indent: int | None = None,
        width: int | None = None,
        allow_unicode: Incomplete | None = None,
        line_break: str | None = None,
        encoding: str | None = None,
        explicit_start: Incomplete | None = None,
        explicit_end: Incomplete | None = None,
        version: Sequence[int] | None = None,
        tags: Mapping[str, str] | None = None,
        sort_keys: bool = True,
    ) -> None: ...

class CDumper(CEmitter, SafeRepresenter, Resolver):
    def __init__(
        self,
        stream: IO[Any],
        default_style: str | None = None,
        default_flow_style: bool = False,
        canonical: Incomplete | None = None,
        indent: int | None = None,
        width: int | None = None,
        allow_unicode: Incomplete | None = None,
        line_break: str | None = None,
        encoding: str | None = None,
        explicit_start: Incomplete | None = None,
        explicit_end: Incomplete | None = None,
        version: Sequence[int] | None = None,
        tags: Mapping[str, str] | None = None,
        sort_keys: bool = True,
    ) -> None: ...

CSafeDumper = CDumper
