from datetime import timedelta
from typing import (
    ClassVar,
    Literal,
    Self,
    TypeAlias,
    overload,
)

import numpy as np

from pandas._libs.tslibs import (
    NaTType,
    Tick,
)
from pandas._typing import (
    Frequency,
    TimeUnit,
    npt,
)

# This should be kept consistent with the keys in the dict timedelta_abbrevs
# in pandas/_libs/tslibs/timedeltas.pyx
UnitChoices: TypeAlias = Literal[
    "Y",
    "y",
    "M",
    "W",
    "w",
    "D",
    "d",
    "days",
    "day",
    "hours",
    "hour",
    "hr",
    "h",
    "m",
    "minute",
    "min",
    "minutes",
    "s",
    "seconds",
    "sec",
    "second",
    "ms",
    "milliseconds",
    "millisecond",
    "milli",
    "millis",
    "us",
    "microseconds",
    "microsecond",
    "µs",
    "micro",
    "micros",
    "ns",
    "nanoseconds",
    "nano",
    "nanos",
    "nanosecond",
]

def get_unit_for_round(freq, creso: int) -> int: ...
def disallow_ambiguous_unit(unit: str | None) -> None: ...
def ints_to_pytimedelta(
    m8values: npt.NDArray[np.timedelta64],
    box: bool = ...,
) -> npt.NDArray[np.object_]: ...
def array_to_timedelta64(
    values: npt.NDArray[np.object_],
    unit: str | None = ...,
    errors: str = ...,
    creso: int = ...,
) -> np.ndarray: ...  # np.ndarray[m8ns]
def parse_timedelta_unit(unit: str | None) -> UnitChoices: ...
def delta_to_nanoseconds(
    delta: np.timedelta64 | timedelta | Tick,
    reso: int = ...,  # NPY_DATETIMEUNIT
    round_ok: bool = ...,
) -> int: ...
def floordiv_object_array(
    left: np.ndarray, right: npt.NDArray[np.object_]
) -> np.ndarray: ...
def truediv_object_array(
    left: np.ndarray, right: npt.NDArray[np.object_]
) -> np.ndarray: ...

# error: Definition of "__eq__" in base class "timedelta" is incompatible with
# definition in base class "NaTType"
class Timedelta(timedelta, NaTType):  # type: ignore[misc]
    _creso: int
    min: ClassVar[Timedelta]
    max: ClassVar[Timedelta]
    resolution: ClassVar[Timedelta]
    value: int  # np.int64
    _value: int  # np.int64
    def __new__(
        cls: type[Self],
        value=...,
        unit: str | None = ...,
        **kwargs: float | np.integer | np.floating,
    ) -> Self: ...
    @classmethod
    def _from_value_and_reso(cls, value: np.int64, reso: int) -> Timedelta: ...
    @property
    def days(self) -> int: ...
    @property
    def seconds(self) -> int: ...
    @property
    def microseconds(self) -> int: ...
    def total_seconds(self) -> float: ...
    def to_pytimedelta(self) -> timedelta: ...
    def to_timedelta64(self) -> np.timedelta64: ...
    # error: Signature of "asm8" incompatible with supertype "NaTType"
    @property
    def asm8(self) -> np.timedelta64: ...  # type: ignore[override]
    # error: Signature of "round" incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"
    def round(self, freq: Frequency | timedelta) -> Self: ...  # type: ignore[override]
    # error: Signature of "floor" incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"
    def floor(self, freq: Frequency | timedelta) -> Self: ...  # type: ignore[override]
    # error: Signature of "ceil" incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"
    def ceil(self, freq: Frequency | timedelta) -> Self: ...  # type: ignore[override]
    @property
    def resolution_string(self) -> str: ...
    # error: Argument 1 of "__add__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument type as
    # "Timedelta | datetime | timedelta | Period | datetime64[date | int | None] |
    # timedelta64[timedelta | int | None]"
    def __add__(self, other: timedelta) -> Timedelta: ...  # type: ignore[override]
    # error: Argument 1 of "__radd__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument type as
    # "Timedelta | datetime | timedelta | Period | datetime64[date | int | None] |
    # timedelta64[timedelta | int | None]"
    def __radd__(self, other: timedelta) -> Timedelta: ...  # type: ignore[override]
    # error: Argument 1 of "__sub__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument type as
    # "Timedelta | datetime | timedelta | Period | datetime64[date | int | None] |
    # timedelta64[timedelta | int | None]"
    def __sub__(self, other: timedelta) -> Timedelta: ...  # type: ignore[override]
    # error: Argument 1 of "__rsub__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument type as
    # "Timedelta | datetime | timedelta | Period | datetime64[date | int | None] |
    # timedelta64[timedelta | int | None]"
    def __rsub__(self, other: timedelta) -> Timedelta: ...  # type: ignore[override]
    def __neg__(self) -> Timedelta: ...
    def __pos__(self) -> Timedelta: ...
    def __abs__(self) -> Timedelta: ...
    def __mul__(self, other: float) -> Timedelta: ...
    def __rmul__(self, other: float) -> Timedelta: ...
    # error: Signature of "__floordiv__" incompatible with supertype "timedelta"
    @overload  # type: ignore[override]
    def __floordiv__(self, other: timedelta) -> int: ...
    @overload
    def __floordiv__(self, other: float) -> Timedelta: ...
    @overload
    def __floordiv__(
        self, other: npt.NDArray[np.timedelta64]
    ) -> npt.NDArray[np.intp]: ...
    @overload
    def __floordiv__(
        self, other: npt.NDArray[np.number]
    ) -> npt.NDArray[np.timedelta64] | Timedelta: ...
    @overload
    def __rfloordiv__(self, other: timedelta | str) -> int: ...
    @overload
    def __rfloordiv__(self, other: None | NaTType) -> NaTType: ...
    @overload
    def __rfloordiv__(self, other: np.ndarray) -> npt.NDArray[np.timedelta64]: ...
    # error: Argument 1 of "__truediv__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument type as
    # "Timedelta | datetime | timedelta | Period | datetime64[date | int | None] |
    # timedelta64[timedelta | int | None]"
    @overload  # type: ignore[override]
    def __truediv__(self, other: Self | timedelta | np.timedelta64, /) -> float: ...
    @overload
    def __truediv__(self, other: float) -> Self: ...
    def __mod__(self, other: timedelta) -> Timedelta: ...
    def __divmod__(self, other: timedelta) -> tuple[int, Timedelta]: ...
    # error: Return type "bool" of "__le__" incompatible with return type
    # "Literal[False]" in supertype "pandas._libs.tslibs.nattype.NaTType"
    # error: Argument 1 of "__le__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument
    # type as "Timedelta | datetime | timedelta | Period | datetime64[date | int
    # | None] | timedelta64[timedelta | int | None]"
    def __le__(self, other: timedelta) -> bool: ...  # type: ignore[override]
    # error: Return type "bool" of "__lt__" incompatible with return type
    # "Literal[False]" in supertype "pandas._libs.tslibs.nattype.NaTType"
    # error: Argument 1 of "__lt__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument
    # type as "Timedelta | datetime | timedelta | Period | datetime64[date | int
    # | None] | timedelta64[timedelta | int | None]"
    def __lt__(self, other: timedelta) -> bool: ...  # type: ignore[override]
    # error: Return type "bool" of "__ge__" incompatible with return type
    # "Literal[False]" in supertype "pandas._libs.tslibs.nattype.NaTType"
    # error: Argument 1 of "__ge__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument
    # type as "Timedelta | datetime | timedelta | Period | datetime64[date | int
    # | None] | timedelta64[timedelta | int | None]"
    def __ge__(self, other: timedelta) -> bool: ...  # type: ignore[override]
    # error: Return type "bool" of "__gt__" incompatible with return type
    # "Literal[False]" in supertype "pandas._libs.tslibs.nattype.NaTType"
    # error: Argument 1 of "__gt__" is incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"; supertype defines the argument
    # type as "Timedelta | datetime | timedelta | Period | datetime64[date | int
    # | None] | timedelta64[timedelta | int | None]"
    def __gt__(self, other: timedelta) -> bool: ...  # type: ignore[override]
    def __hash__(self) -> int: ...
    #  error: Signature of "isoformat" incompatible with supertype
    # "pandas._libs.tslibs.nattype.NaTType"
    def isoformat(self) -> str: ...  # type: ignore[override]
    # error: Return type "timedelta64[timedelta | int | None]" of "to_numpy"
    # incompatible with return type "datetime64[date | int | None]" in supertype
    # "pandas._libs.tslibs.nattype.NaTType"
    def to_numpy(  # type: ignore[override]
        self, dtype: npt.DTypeLike | None = ..., copy: bool = False
    ) -> np.timedelta64: ...
    def view(self, dtype: npt.DTypeLike) -> object: ...
    @property
    def unit(self) -> TimeUnit: ...
    def as_unit(self, unit: TimeUnit, round_ok: bool = ...) -> Timedelta: ...
