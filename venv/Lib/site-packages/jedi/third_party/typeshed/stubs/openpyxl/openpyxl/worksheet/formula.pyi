from _typeshed import Incomplete
from collections.abc import Iterator
from typing import ClassVar

class DataTableFormula:
    t: ClassVar[str]

    ref: Incomplete
    ca: bool
    dt2D: bool
    dtr: bool
    r1: Incomplete | None
    r2: Incomplete | None
    del1: bool
    del2: bool

    def __init__(
        self,
        ref,
        ca: bool = False,
        dt2D: bool = False,
        dtr: bool = False,
        r1=None,
        r2=None,
        del1: bool = False,
        del2: bool = False,
        **kw,
    ) -> None: ...
    def __iter__(self) -> Iterator[tuple[str, str]]: ...

class ArrayFormula:
    t: ClassVar[str]
    ref: Incomplete
    text: Incomplete | None

    def __init__(self, ref, text=None) -> None: ...
    def __iter__(self) -> Iterator[tuple[str, str]]: ...
