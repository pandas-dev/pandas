from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, ReadableBuffer, Unused
from collections.abc import Iterable, Sized
from datetime import datetime
from re import Pattern
from typing import Any, Generic, Literal, TypeVar, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors import Strict
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.fill import Blip
from openpyxl.worksheet.cell_range import CellRange, MultiCellRange

_T = TypeVar("_T")
_P = TypeVar("_P", str, ReadableBuffer)
_N = TypeVar("_N", bound=bool, default=Literal[False])
_L = TypeVar("_L", bound=Sized)
_M = TypeVar("_M", int, float)

_ExpectedTypeParam: TypeAlias = type[_T] | tuple[type[_T], ...]
_ConvertibleToMultiCellRange: TypeAlias = MultiCellRange | str | Iterable[CellRange]
# Since everything is convertible to a bool, this restricts to only intended expected types of intended literals
_ConvertibleToBool: TypeAlias = bool | str | int | None  # True | False | "true" | "t" | "false" | "f" | 1 | 0 | None

class Descriptor(Generic[_T]):
    name: str | None
    def __init__(self, name: str | None = None, **kw: object) -> None: ...
    def __get__(self, instance: Serialisable | Strict, cls: type | None) -> _T: ...
    def __set__(self, instance: Serialisable | Strict, value: _T) -> None: ...

class Typed(Descriptor[_T], Generic[_T, _N]):
    __doc__: str
    # Members optional in __init__
    expected_type: type[_T]
    allow_none: _N
    nested: bool

    @overload
    def __init__(
        self: Typed[_T, Literal[True]],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: Literal[True],
        nested: bool = False,
    ) -> None: ...
    @overload
    def __init__(
        self: Typed[_T, Literal[False]],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: Literal[False] = False,
        nested: bool = False,
    ) -> None: ...
    @overload
    def __get__(self: Typed[_T, Literal[True]], instance: Serialisable | Strict, cls: type | None = None) -> _T | None: ...
    @overload
    def __get__(self: Typed[_T, Literal[False]], instance: Serialisable | Strict, cls: type | None = None) -> _T: ...
    @overload
    def __set__(self: Typed[_T, Literal[True]], instance: Serialisable | Strict, value: _T | None) -> None: ...
    @overload
    def __set__(self: Typed[_T, Literal[False]], instance: Serialisable | Strict, value: _T) -> None: ...

class Convertible(Typed[_T, _N]):
    @overload
    def __init__(
        self: Convertible[_T, Literal[True]],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: Literal[True],
    ) -> None: ...
    @overload
    def __init__(
        self: Convertible[_T, Literal[False]],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: Literal[False] = False,
    ) -> None: ...
    # NOTE: It is currently impossible to make a generic based on the parameter type of another generic
    # So we implement explicitly the types used internally
    # MultiCellRange
    @overload
    def __set__(
        self: Convertible[MultiCellRange, Literal[True]],
        instance: Serialisable | Strict,
        value: _ConvertibleToMultiCellRange | None,
    ) -> None: ...
    @overload
    def __set__(
        self: Convertible[MultiCellRange, Literal[False]], instance: Serialisable | Strict, value: _ConvertibleToMultiCellRange
    ) -> None: ...
    # str | Blip
    @overload
    def __set__(
        self: Convertible[str, _N] | Convertible[Blip, _N],
        instance: Serialisable | Strict,
        value: object,  # Not[None] when _N = False
    ) -> None: ...
    # bool
    @overload
    def __set__(self: Convertible[bool, _N], instance: Serialisable | Strict, value: _ConvertibleToBool) -> None: ...
    # int
    @overload
    def __set__(
        self: Convertible[int, Literal[True]], instance: Serialisable | Strict, value: ConvertibleToInt | None
    ) -> None: ...
    @overload
    def __set__(self: Convertible[int, Literal[False]], instance: Serialisable | Strict, value: ConvertibleToInt) -> None: ...
    # float
    @overload
    def __set__(
        self: Convertible[float, Literal[True]], instance: Serialisable | Strict, value: ConvertibleToFloat | None
    ) -> None: ...
    @overload
    def __set__(self: Convertible[float, Literal[False]], instance: Serialisable | Strict, value: ConvertibleToFloat) -> None: ...
    # Anything else
    @overload
    def __set__(self: Convertible[_T, Literal[True]], instance: Serialisable | Strict, value: _T | int | Any | None) -> None: ...

class Max(Convertible[_M, _N]):
    expected_type: type[_M]
    allow_none: _N
    max: float
    @overload
    def __init__(
        self: Max[int, Literal[True]], *, expected_type: _ExpectedTypeParam[int], allow_none: Literal[True], max: float
    ) -> None: ...
    @overload
    def __init__(
        self: Max[int, Literal[False]], *, expected_type: _ExpectedTypeParam[int], allow_none: Literal[False] = False, max: float
    ) -> None: ...
    # mypy can't infer type from `expected_type = float` (pyright can), so we have to add extra overloads
    @overload
    def __init__(
        self: Max[float, Literal[True]], *, expected_type: _ExpectedTypeParam[float] = ..., allow_none: Literal[True], max: float
    ) -> None: ...
    @overload
    def __init__(
        self: Max[float, Literal[False]],
        *,
        expected_type: _ExpectedTypeParam[float] = ...,
        allow_none: Literal[False] = False,
        max: float,
    ) -> None: ...
    @overload  # type: ignore[override]  # Different restrictions
    def __set__(self: Max[int, Literal[True]], instance: Serialisable | Strict, value: ConvertibleToInt | None) -> None: ...
    @overload
    def __set__(self: Max[int, Literal[False]], instance: Serialisable | Strict, value: ConvertibleToInt) -> None: ...
    @overload
    def __set__(self: Max[float, Literal[True]], instance: Serialisable | Strict, value: ConvertibleToFloat | None) -> None: ...
    @overload
    def __set__(self: Max[float, Literal[False]], instance: Serialisable | Strict, value: ConvertibleToFloat) -> None: ...

class Min(Convertible[_M, _N]):
    expected_type: type[_M]
    allow_none: _N
    min: float
    @overload
    def __init__(
        self: Min[int, Literal[True]], *, expected_type: _ExpectedTypeParam[int], allow_none: Literal[True], min: float
    ) -> None: ...
    @overload
    def __init__(
        self: Min[int, Literal[False]], *, expected_type: _ExpectedTypeParam[int], allow_none: Literal[False] = False, min: float
    ) -> None: ...
    # mypy can't infer type from `expected_type = float` (pyright can), so we have to add extra overloads
    @overload
    def __init__(
        self: Min[float, Literal[True]], *, expected_type: _ExpectedTypeParam[float] = ..., allow_none: Literal[True], min: float
    ) -> None: ...
    @overload
    def __init__(
        self: Min[float, Literal[False]],
        *,
        expected_type: _ExpectedTypeParam[float] = ...,
        allow_none: Literal[False] = False,
        min: float,
    ) -> None: ...
    @overload  # type: ignore[override]  # Different restrictions
    def __set__(self: Min[int, Literal[True]], instance: Serialisable | Strict, value: ConvertibleToInt | None) -> None: ...
    @overload
    def __set__(self: Min[int, Literal[False]], instance: Serialisable | Strict, value: ConvertibleToInt) -> None: ...
    @overload
    def __set__(self: Min[float, Literal[True]], instance: Serialisable | Strict, value: ConvertibleToFloat | None) -> None: ...
    @overload
    def __set__(self: Min[float, Literal[False]], instance: Serialisable | Strict, value: ConvertibleToFloat) -> None: ...

class MinMax(Min[_M, _N], Max[_M, _N]):
    expected_type: type[_M]
    allow_none: _N
    @overload
    def __init__(
        self: MinMax[int, Literal[True]],
        *,
        expected_type: _ExpectedTypeParam[int],
        allow_none: Literal[True],
        min: float,
        max: float,
    ) -> None: ...
    @overload
    def __init__(
        self: MinMax[int, Literal[False]],
        *,
        expected_type: _ExpectedTypeParam[int],
        allow_none: Literal[False] = False,
        min: float,
        max: float,
    ) -> None: ...
    # mypy can't infer type from `expected_type = float` (pyright can), so we have to add extra overloads
    @overload
    def __init__(
        self: MinMax[float, Literal[True]],
        *,
        expected_type: _ExpectedTypeParam[float] = ...,
        allow_none: Literal[True],
        min: float,
        max: float,
    ) -> None: ...
    @overload
    def __init__(
        self: MinMax[float, Literal[False]],
        *,
        expected_type: _ExpectedTypeParam[float] = ...,
        allow_none: Literal[False] = False,
        min: float,
        max: float,
    ) -> None: ...

class Set(Descriptor[_T]):
    __doc__: str
    values: Iterable[_T]
    def __init__(self, name: str | None = None, *, values: Iterable[_T]) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value: _T) -> None: ...

class NoneSet(Set[_T | None]):
    def __init__(self, name: str | None = None, *, values: Iterable[_T | None]) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value: _T | Literal["none"] | None) -> None: ...

class Integer(Convertible[int, _N]):
    allow_none: _N
    expected_type: type[int]
    @overload
    def __init__(self: Integer[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: Integer[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...

class Float(Convertible[float, _N]):
    allow_none: _N
    expected_type: type[float]
    @overload
    def __init__(self: Float[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: Float[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...

class Bool(Convertible[bool, _N]):
    expected_type: type[bool]
    allow_none: _N
    @overload
    def __init__(self: Bool[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: Bool[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value: _ConvertibleToBool) -> None: ...

class String(Typed[str, _N]):
    allow_none: _N
    expected_type: type[str]
    @overload
    def __init__(self: String[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: String[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...

class Text(String[_N], Convertible[str, _N]): ...  # unused

class ASCII(Typed[bytes, _N]):  # unused
    expected_type: type[bytes]
    def __init__(self, name: str | None = None, *, allow_none: bool = False) -> None: ...

class Tuple(Typed[tuple[Any, ...], _N]):  # unused
    expected_type: type[tuple[Any, ...]]
    def __init__(self, name: str | None = None, *, allow_none: bool = False) -> None: ...

class Length(Descriptor[_L]):
    def __init__(self, name: Unused = None, *, length: int) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value: _L) -> None: ...

class Default(Typed[_T, _N]):  # unused
    def __init__(
        self, name: Unused = None, *, expected_type: _ExpectedTypeParam[_T], allow_none: bool = False, defaults: Unused = {}
    ) -> None: ...
    def __call__(self) -> _T: ...

# Note: Aliases types can't be inferred. Anyway an alias means there's another option.
# Incomplete: Make it generic with explicit getter/setter type arguments?
class Alias(Descriptor[Incomplete]):
    alias: str
    def __init__(self, alias: str) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value) -> None: ...
    def __get__(self, instance: Serialisable | Strict, cls: Unused): ...

class MatchPattern(Descriptor[_P], Generic[_P, _N]):
    allow_none: _N
    test_pattern: Pattern[bytes] | Pattern[str]
    pattern: str | Pattern[str] | bytes | Pattern[bytes]

    @overload  # str
    def __init__(
        self: MatchPattern[str, Literal[True]], name: str | None = None, *, pattern: str | Pattern[str], allow_none: Literal[True]
    ) -> None: ...
    @overload  # str | None
    def __init__(
        self: MatchPattern[str, Literal[False]],
        name: str | None = None,
        *,
        pattern: str | Pattern[str],
        allow_none: Literal[False] = False,
    ) -> None: ...
    @overload  # bytes
    def __init__(
        self: MatchPattern[ReadableBuffer, Literal[True]],
        name: str | None = None,
        *,
        pattern: bytes | Pattern[bytes],
        allow_none: Literal[True],
    ) -> None: ...
    @overload  # bytes | None
    def __init__(
        self: MatchPattern[ReadableBuffer, Literal[False]],
        name: str | None = None,
        *,
        pattern: bytes | Pattern[bytes],
        allow_none: Literal[False] = False,
    ) -> None: ...
    @overload
    def __get__(self: MatchPattern[_P, Literal[True]], instance: Serialisable | Strict, cls: type | None = None) -> _P | None: ...
    @overload
    def __get__(self: MatchPattern[_P, Literal[False]], instance: Serialisable | Strict, cls: type | None = None) -> _P: ...
    @overload
    def __set__(self: MatchPattern[_P, Literal[True]], instance: Serialisable | Strict, value: _P | None) -> None: ...
    @overload
    def __set__(self: MatchPattern[_P, Literal[False]], instance: Serialisable | Strict, value: _P) -> None: ...

class DateTime(Typed[datetime, _N]):
    allow_none: _N
    expected_type: type[datetime]
    @overload
    def __init__(self: DateTime[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: DateTime[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...
    @overload
    def __set__(self: DateTime[Literal[True]], instance: Serialisable | Strict, value: datetime | str | None) -> None: ...
    @overload
    def __set__(self: DateTime[Literal[False]], instance: Serialisable | Strict, value: datetime | str) -> None: ...
