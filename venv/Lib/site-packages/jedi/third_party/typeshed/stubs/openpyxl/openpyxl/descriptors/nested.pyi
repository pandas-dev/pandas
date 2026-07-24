from _typeshed import ConvertibleToFloat, ConvertibleToInt, Unused
from collections.abc import Iterable
from typing import Any, ClassVar, Literal, NoReturn, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors import Strict
from openpyxl.descriptors.base import Bool, Convertible, Descriptor, Float, Integer, MinMax, NoneSet, Set, String
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.fill import Blip
from openpyxl.xml.functions import Element

from ..xml._functions_overloads import _HasGet, _HasTagAndGet, _HasText
from .base import _M, _N, _T, _ConvertibleToBool, _ExpectedTypeParam

_NestedNoneSetParam: TypeAlias = _HasTagAndGet[_T | Literal["none"] | None] | _T | Literal["none"] | None

# NOTE: type: ignore[misc]: Class does not reimplement the relevant methods, so runtime also has incompatible supertypes

class Nested(Descriptor[_T]):
    nested: ClassVar[Literal[True]]
    attribute: ClassVar[str]
    # Members optional in __init__
    expected_type: type[_T]
    allow_none: bool
    namespace: str | None
    # In usage, "Nested" is closed to "Typed" than "Descriptor", but doesn't use allow_none
    def __init__(
        self: Nested[_T],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: bool = False,
        nested: Unused = True,
        namespace: str | None = None,
    ) -> None: ...
    def __get__(self, instance: Serialisable | Strict, cls: type | None) -> _T: ...
    def __set__(self, instance: Serialisable | Strict, value: _HasTagAndGet[_T] | _T) -> None: ...
    def from_tree(self, node: _HasGet[_T]) -> _T: ...
    @overload
    def to_tree(self, tagname: Unused = None, value: None = None, namespace: Unused = None) -> None: ...
    @overload
    def to_tree(self, tagname: str, value: object, namespace: str | None = None) -> Element: ...

class NestedValue(Nested[_T], Convertible[_T, _N]):  # type: ignore[misc]
    @overload
    def __init__(
        self: NestedValue[_T, Literal[True]],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: Literal[True],
    ) -> None: ...
    @overload
    def __init__(
        self: NestedValue[_T, Literal[False]],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: Literal[False] = False,
    ) -> None: ...
    @overload
    def __get__(self: NestedValue[_T, Literal[True]], instance: Serialisable | Strict, cls: type | None = None) -> _T | None: ...
    @overload
    def __get__(self: NestedValue[_T, Literal[False]], instance: Serialisable | Strict, cls: type | None = None) -> _T: ...
    # NOTE: It is currently impossible to make a generic based on the parameter type of another generic
    # So we implement explicitly the types used internally
    # str | Blip
    @overload
    def __set__(
        self: NestedValue[str, _N] | NestedValue[Blip, _N],
        instance: Serialisable | Strict,
        value: object,  # Not[None] when _N = False
    ) -> None: ...
    # bool
    @overload
    def __set__(
        self: NestedValue[bool, _N],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool,
    ) -> None: ...
    # int
    @overload
    def __set__(
        self: NestedValue[int, Literal[True]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
    ) -> None: ...
    @overload
    def __set__(
        self: NestedValue[int, Literal[False]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt,
    ) -> None: ...
    # float
    @overload
    def __set__(
        self: NestedValue[float, Literal[True]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None,
    ) -> None: ...
    @overload
    def __set__(
        self: NestedValue[float, Literal[False]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[ConvertibleToFloat] | ConvertibleToFloat,
    ) -> None: ...
    # Anything else
    @overload
    def __set__(
        self: NestedValue[_T, Literal[True]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[_T | int | Any] | _T | int | Any | None,
    ) -> None: ...

class NestedText(NestedValue[_T, _N]):
    @overload
    def __init__(
        self: NestedText[_T, Literal[True]],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: Literal[True],
    ) -> None: ...
    @overload
    def __init__(
        self: NestedText[_T, Literal[False]],  # pyright: ignore[reportInvalidTypeVarUse]  #11780
        name: str | None = None,
        *,
        expected_type: _ExpectedTypeParam[_T],
        allow_none: Literal[False] = False,
    ) -> None: ...
    @overload
    def __get__(self: NestedText[_T, Literal[True]], instance: Serialisable | Strict, cls: type | None = None) -> _T | None: ...
    @overload
    def __get__(self: NestedText[_T, Literal[False]], instance: Serialisable | Strict, cls: type | None = None) -> _T: ...
    # NOTE: It is currently impossible to make a generic based on the parameter type of another generic
    # So we implement explicitly the types used internally
    # str
    @overload
    def __set__(  # type: ignore[overload-overlap]
        self: NestedText[str, _N], instance: Serialisable | Strict, value: object  # Not[None] when _N = False
    ) -> None: ...
    # int
    @overload
    def __set__(
        self: NestedText[int, Literal[True]], instance: Serialisable | Strict, value: ConvertibleToInt | None
    ) -> None: ...
    @overload
    def __set__(self: NestedText[int, Literal[False]], instance: Serialisable | Strict, value: ConvertibleToInt) -> None: ...
    # If expected type (_T) is not str, it's impossible to use an Element as the value
    @overload
    def __set__(self: NestedText[_T, Literal[True]], instance: Serialisable | Strict, value: _HasTagAndGet[Any]) -> NoReturn: ...
    # Anything else
    @overload
    def __set__(self: NestedText[_T, Literal[True]], instance: Serialisable | Strict, value: _T | int | Any | None) -> None: ...
    def from_tree(self, node: _HasText) -> str: ...  # type: ignore[override]
    @overload
    def to_tree(self, tagname: Unused = None, value: None = None, namespace: Unused = None) -> None: ...
    @overload
    def to_tree(self, tagname: str, value: object, namespace: str | None = None) -> Element: ...

class NestedFloat(NestedValue[float, _N], Float[_N]):  # type: ignore[misc]
    @overload
    def __init__(self: NestedFloat[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: NestedFloat[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...

class NestedInteger(NestedValue[int, _N], Integer[_N]):  # type: ignore[misc]
    @overload
    def __init__(self: NestedInteger[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: NestedInteger[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...

class NestedString(NestedValue[str, _N], String[_N]):  # type: ignore[misc]
    @overload
    def __init__(self: NestedString[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: NestedString[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...

class NestedBool(NestedValue[bool, _N], Bool[_N]):  # type: ignore[misc]
    @overload
    def __init__(self: NestedBool[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: NestedBool[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool) -> None: ...
    def from_tree(self, node: _HasGet[bool]) -> bool: ...

class NestedNoneSet(Nested[_T | None], NoneSet[_T]):
    def __init__(self, name: str | None = None, *, values: Iterable[_T | None]) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value: _NestedNoneSetParam[_T]) -> None: ...

class NestedSet(Nested[_T], Set[_T]):
    def __init__(self, name: str | None = None, *, values: Iterable[_T]) -> None: ...

class NestedMinMax(Nested[_M], MinMax[_M, _N]):  # type: ignore[misc]
    @overload
    def __init__(
        self: NestedMinMax[int, Literal[True]],
        *,
        expected_type: _ExpectedTypeParam[int],
        allow_none: Literal[True],
        min: float,
        max: float,
    ) -> None: ...
    @overload
    def __init__(
        self: NestedMinMax[int, Literal[False]],
        *,
        expected_type: _ExpectedTypeParam[int],
        allow_none: Literal[False] = False,
        min: float,
        max: float,
    ) -> None: ...
    # mypy can't infer type from `expected_type = float` (pyright can), so we have to add extra overloads
    @overload
    def __init__(
        self: NestedMinMax[float, Literal[True]],
        *,
        expected_type: _ExpectedTypeParam[float] = ...,
        allow_none: Literal[True],
        min: float,
        max: float,
    ) -> None: ...
    @overload
    def __init__(
        self: NestedMinMax[float, Literal[False]],
        *,
        expected_type: _ExpectedTypeParam[float] = ...,
        allow_none: Literal[False] = False,
        min: float,
        max: float,
    ) -> None: ...
    @overload
    def __get__(self: NestedMinMax[_M, Literal[True]], instance: Serialisable | Strict, cls: type | None = None) -> _M | None: ...
    @overload
    def __get__(self: NestedMinMax[_M, Literal[False]], instance: Serialisable | Strict, cls: type | None = None) -> _M: ...
    @overload  # type: ignore[override]  # Different restrictions
    def __set__(
        self: NestedMinMax[int, Literal[True]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
    ) -> None: ...
    @overload
    def __set__(
        self: NestedMinMax[int, Literal[False]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt,
    ) -> None: ...
    @overload
    def __set__(
        self: NestedMinMax[float, Literal[True]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None,
    ) -> None: ...
    @overload
    def __set__(
        self: NestedMinMax[float, Literal[False]],
        instance: Serialisable | Strict,
        value: _HasTagAndGet[ConvertibleToFloat] | ConvertibleToFloat,
    ) -> None: ...

class EmptyTag(Nested[bool], Bool[_N]):  # type: ignore[misc]
    @overload
    def __init__(self: EmptyTag[Literal[True]], name: str | None = None, *, allow_none: Literal[True]) -> None: ...
    @overload
    def __init__(self: EmptyTag[Literal[False]], name: str | None = None, *, allow_none: Literal[False] = False) -> None: ...
    @overload
    def __get__(self: EmptyTag[Literal[True]], instance: Serialisable | Strict, cls: type | None = None) -> bool | None: ...
    @overload
    def __get__(self: EmptyTag[Literal[False]], instance: Serialisable | Strict, cls: type | None = None) -> bool: ...
    def __set__(self, instance: Serialisable | Strict, value: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool) -> None: ...
    def from_tree(self, node: Unused) -> Literal[True]: ...
    @overload
    def to_tree(self, tagname: Unused = None, value: None = None, namespace: Unused = None) -> None: ...
    @overload
    def to_tree(self, tagname: str, value: object, namespace: str | None = None) -> Element: ...
