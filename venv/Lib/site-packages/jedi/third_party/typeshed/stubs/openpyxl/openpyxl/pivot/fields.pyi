from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete
from datetime import datetime
from typing import ClassVar, Literal, overload

from openpyxl.descriptors.base import Bool, DateTime, Float, Integer, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.serialisable import Serialisable

class Index(Serialisable):
    tagname: ClassVar[str]
    v: Integer[Literal[True]]
    def __init__(self, v: ConvertibleToInt | None = 0) -> None: ...

class Tuple(Serialisable):
    fld: Integer[Literal[True]]
    hier: Integer[Literal[True]]
    item: Integer[Literal[False]]
    def __init__(self, fld: ConvertibleToInt, hier: ConvertibleToInt, item: ConvertibleToInt) -> None: ...

class TupleList(Serialisable):
    c: Integer[Literal[True]]
    tpl: Typed[Tuple, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(self, c: ConvertibleToInt | None = None, *, tpl: Tuple) -> None: ...
    @overload
    def __init__(self, c: ConvertibleToInt | None, tpl: Tuple) -> None: ...

class Missing(Serialisable):
    tagname: ClassVar[str]
    tpls: Incomplete
    x: Incomplete
    u: Bool[Literal[True]]
    f: Bool[Literal[True]]
    c: String[Literal[True]]
    cp: Integer[Literal[True]]
    _in: Integer[Literal[True]]  # Not private. Avoids name clash
    bc: Incomplete
    fc: Incomplete
    i: Bool[Literal[True]]
    un: Bool[Literal[True]]
    st: Bool[Literal[True]]
    b: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        tpls=(),
        x=(),
        u: _ConvertibleToBool | None = None,
        f: _ConvertibleToBool | None = None,
        c: str | None = None,
        cp: ConvertibleToInt | None = None,
        _in: ConvertibleToInt | None = None,
        bc=None,
        fc=None,
        i: _ConvertibleToBool | None = None,
        un: _ConvertibleToBool | None = None,
        st: _ConvertibleToBool | None = None,
        b: _ConvertibleToBool | None = None,
    ) -> None: ...

class Number(Serialisable):
    tagname: ClassVar[str]
    tpls: Incomplete
    x: Incomplete
    v: Float[Literal[False]]
    u: Bool[Literal[True]]
    f: Bool[Literal[True]]
    c: String[Literal[True]]
    cp: Integer[Literal[True]]
    _in: Integer[Literal[True]]  # Not private. Avoids name clash
    bc: Incomplete
    fc: Incomplete
    i: Bool[Literal[True]]
    un: Bool[Literal[True]]
    st: Bool[Literal[True]]
    b: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        tpls=(),
        x=(),
        *,
        v: ConvertibleToFloat,
        u: _ConvertibleToBool | None = None,
        f: _ConvertibleToBool | None = None,
        c: str | None = None,
        cp: ConvertibleToInt | None = None,
        _in: ConvertibleToInt | None = None,
        bc=None,
        fc=None,
        i: _ConvertibleToBool | None = None,
        un: _ConvertibleToBool | None = None,
        st: _ConvertibleToBool | None = None,
        b: _ConvertibleToBool | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        tpls,
        x,
        v: ConvertibleToFloat,
        u: _ConvertibleToBool | None = None,
        f: _ConvertibleToBool | None = None,
        c: str | None = None,
        cp: ConvertibleToInt | None = None,
        _in: ConvertibleToInt | None = None,
        bc=None,
        fc=None,
        i: _ConvertibleToBool | None = None,
        un: _ConvertibleToBool | None = None,
        st: _ConvertibleToBool | None = None,
        b: _ConvertibleToBool | None = None,
    ) -> None: ...

class Error(Serialisable):
    tagname: ClassVar[str]
    tpls: Typed[TupleList, Literal[True]]
    x: Incomplete
    v: String[Literal[False]]
    u: Bool[Literal[True]]
    f: Bool[Literal[True]]
    c: String[Literal[True]]
    cp: Integer[Literal[True]]
    _in: Integer[Literal[True]]  # Not private. Avoids name clash
    bc: Incomplete
    fc: Incomplete
    i: Bool[Literal[True]]
    un: Bool[Literal[True]]
    st: Bool[Literal[True]]
    b: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        tpls: TupleList | None = None,
        x=(),
        *,
        v: str,
        u: _ConvertibleToBool | None = None,
        f: _ConvertibleToBool | None = None,
        c: str | None = None,
        cp: ConvertibleToInt | None = None,
        _in: ConvertibleToInt | None = None,
        bc=None,
        fc=None,
        i: _ConvertibleToBool | None = None,
        un: _ConvertibleToBool | None = None,
        st: _ConvertibleToBool | None = None,
        b: _ConvertibleToBool | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        tpls: TupleList | None,
        x,
        v: str,
        u: _ConvertibleToBool | None = None,
        f: _ConvertibleToBool | None = None,
        c: str | None = None,
        cp: ConvertibleToInt | None = None,
        _in: ConvertibleToInt | None = None,
        bc=None,
        fc=None,
        i: _ConvertibleToBool | None = None,
        un: _ConvertibleToBool | None = None,
        st: _ConvertibleToBool | None = None,
        b: _ConvertibleToBool | None = None,
    ) -> None: ...

class Boolean(Serialisable):
    tagname: ClassVar[str]
    x: Incomplete
    v: Bool[Literal[False]]
    u: Bool[Literal[True]]
    f: Bool[Literal[True]]
    c: String[Literal[True]]
    cp: Integer[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        x=(),
        v: _ConvertibleToBool = None,
        u: _ConvertibleToBool | None = None,
        f: _ConvertibleToBool | None = None,
        c: str | None = None,
        cp: ConvertibleToInt | None = None,
    ) -> None: ...

class Text(Serialisable):
    tagname: ClassVar[str]
    tpls: Incomplete
    x: Incomplete
    v: String[Literal[False]]
    u: Bool[Literal[True]]
    f: Bool[Literal[True]]
    c: String[Literal[True]]
    cp: Integer[Literal[True]]
    _in: Integer[Literal[True]]  # Not private. Avoids name clash
    bc: Incomplete
    fc: Incomplete
    i: Bool[Literal[True]]
    un: Bool[Literal[True]]
    st: Bool[Literal[True]]
    b: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        tpls=(),
        x=(),
        v=None,
        u: _ConvertibleToBool | None = None,
        f: _ConvertibleToBool | None = None,
        c=None,
        cp: ConvertibleToInt | None = None,
        _in: ConvertibleToInt | None = None,
        bc=None,
        fc=None,
        i: _ConvertibleToBool | None = None,
        un: _ConvertibleToBool | None = None,
        st: _ConvertibleToBool | None = None,
        b: _ConvertibleToBool | None = None,
    ) -> None: ...

class DateTimeField(Serialisable):
    tagname: ClassVar[str]
    x: Incomplete
    v: DateTime[Literal[False]]
    u: Bool[Literal[True]]
    f: Bool[Literal[True]]
    c: String[Literal[True]]
    cp: Integer[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        x=(),
        v: datetime | str | None = None,
        u: _ConvertibleToBool | None = None,
        f: _ConvertibleToBool | None = None,
        c: str | None = None,
        cp: ConvertibleToInt | None = None,
    ) -> None: ...
