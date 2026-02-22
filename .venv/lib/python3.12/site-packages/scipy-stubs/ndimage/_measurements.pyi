from _typeshed import Incomplete
from collections.abc import Callable, Sequence
from typing import Any, Protocol, TypeAlias, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc
from numpy._typing import _ArrayLike, _DTypeLike

__all__ = [
    "center_of_mass",
    "extrema",
    "find_objects",
    "histogram",
    "label",
    "labeled_comprehension",
    "maximum",
    "maximum_position",
    "mean",
    "median",
    "minimum",
    "minimum_position",
    "standard_deviation",
    "sum",
    "sum_labels",
    "value_indices",
    "variance",
    "watershed_ift",
]

_SCT = TypeVar("_SCT", bound=npc.number | np.bool_)
_ISCT = TypeVar("_ISCT", bound=npc.inexact, default=npc.inexact)

__Func1: TypeAlias = Callable[[onp.ToComplex | onp.ToComplexND], onp.ToComplex]
__Func2: TypeAlias = Callable[[onp.ToComplex | onp.ToComplexND, onp.ToComplex | onp.ToComplexND], onp.ToComplex]
_ComprehensionFunc: TypeAlias = __Func1 | __Func2

_Idx0D: TypeAlias = tuple[np.intp, ...]
_IdxND: TypeAlias = list[_Idx0D]

_Extrema0D: TypeAlias = tuple[_SCT, _SCT, _Idx0D, _Idx0D]
_ExtremaND: TypeAlias = tuple[onp.ArrayND[_SCT], onp.ArrayND[_SCT], _IdxND, _IdxND]

_Coord0D: TypeAlias = tuple[np.float64, ...]
_Coord1D: TypeAlias = list[_Coord0D]
_CoordND: TypeAlias = list[tuple[onp.ArrayND[np.float64], ...]]

#
def label(
    input: onp.ToComplex | onp.ToComplexND,
    structure: onp.ToComplex | onp.ToComplexND | None = None,
    output: onp.ArrayND[np.int32 | np.intp] | None = None,
) -> int | tuple[onp.ArrayND[np.int32 | np.intp], int]: ...

#
def find_objects(input: onp.ToInt | onp.ToIntND, max_label: int = 0) -> list[tuple[slice, ...]]: ...

#
def value_indices(
    arr: onp.ToInt | onp.ToIntND, *, ignore_value: int | None = None
) -> dict[np.intp, tuple[onp.ArrayND[np.intp], ...]]: ...

#
@overload
def labeled_comprehension(
    input: onp.ToComplex | onp.ToComplexND,
    labels: onp.ToComplex | onp.ToComplexND | None,
    index: onp.ToInt | onp.ToIntND | None,
    func: _ComprehensionFunc,
    out_dtype: _DTypeLike[_SCT],
    default: onp.ToFloat,
    pass_positions: bool = False,
) -> onp.ArrayND[_SCT]: ...
@overload
def labeled_comprehension(
    input: onp.ToComplex | onp.ToComplexND,
    labels: onp.ToComplex | onp.ToComplexND | None,
    index: onp.ToInt | onp.ToIntND | None,
    func: _ComprehensionFunc,
    out_dtype: onp.AnyIntPDType,
    default: onp.ToInt,
    pass_positions: bool = False,
) -> onp.ArrayND[np.intp]: ...
@overload
def labeled_comprehension(
    input: onp.ToComplex | onp.ToComplexND,
    labels: onp.ToComplex | onp.ToComplexND | None,
    index: onp.ToInt | onp.ToIntND | None,
    func: _ComprehensionFunc,
    out_dtype: onp.AnyFloat64DType | None,
    default: onp.ToFloat,
    pass_positions: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def labeled_comprehension(
    input: onp.ToComplex | onp.ToComplexND,
    labels: onp.ToComplex | onp.ToComplexND | None,
    index: onp.ToInt | onp.ToIntND | None,
    func: _ComprehensionFunc,
    out_dtype: onp.AnyComplex128DType,
    default: onp.ToComplex,
    pass_positions: bool = False,
) -> onp.ArrayND[np.complex128]: ...

#
@type_check_only
class _DefStatistic(Protocol):
    @overload
    def __call__(
        self, /, input: onp.CanArrayND[_ISCT], labels: onp.ToInt | onp.ToIntND | None = None, index: None = None
    ) -> _ISCT: ...
    @overload
    def __call__(
        self, /, input: onp.CanArrayND[_ISCT], labels: onp.ToInt | onp.ToIntND, index: onp.ArrayND[npc.integer]
    ) -> onp.ArrayND[_ISCT]: ...
    @overload
    def __call__(
        self, /, input: onp.CanArrayND[_ISCT], labels: onp.ToInt | onp.ToIntND | None = None, *, index: onp.ArrayND[npc.integer]
    ) -> onp.ArrayND[_ISCT]: ...
    @overload
    def __call__(
        self, /, input: onp.CanArrayND[_ISCT], labels: onp.ToInt | onp.ToIntND, index: onp.ToInt | onp.ToIntND
    ) -> onp.ArrayND[_ISCT]: ...
    @overload
    def __call__(
        self, /, input: onp.CanArrayND[_ISCT], labels: onp.ToInt | onp.ToIntND | None = None, *, index: onp.ToInt | onp.ToIntND
    ) -> onp.ArrayND[_ISCT]: ...
    @overload
    def __call__(
        self, /, input: onp.ToInt | onp.ToIntND, labels: onp.ToInt | onp.ToIntND | None = None, index: None = None
    ) -> np.float64: ...
    @overload
    def __call__(
        self, /, input: onp.ToInt | onp.ToIntND, labels: onp.ToInt | onp.ToIntND, index: onp.ArrayND[npc.integer]
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def __call__(
        self, /, input: onp.ToInt | onp.ToIntND, labels: onp.ToInt | onp.ToIntND | None = None, *, index: onp.ArrayND[npc.integer]
    ) -> onp.ArrayND[np.float64]: ...
    @overload
    def __call__(
        self, /, input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND | None = None, index: None = None
    ) -> npc.inexact: ...
    @overload
    def __call__(
        self, /, input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND, index: onp.ArrayND[npc.integer]
    ) -> onp.ArrayND[npc.inexact]: ...
    @overload
    def __call__(
        self,
        /,
        input: onp.ToComplex | onp.ToComplexND,
        labels: onp.ToInt | onp.ToIntND | None = None,
        *,
        index: onp.ArrayND[npc.integer],
    ) -> onp.ArrayND[npc.inexact]: ...
    @overload
    def __call__(
        self,
        /,
        input: onp.ToComplex | onp.ToComplexND,
        labels: onp.ToInt | onp.ToIntND | None = None,
        index: onp.ToInt | onp.ToIntND | None = None,
    ) -> npc.inexact | onp.ArrayND[npc.inexact]: ...

sum: _DefStatistic
sum_labels: _DefStatistic
mean: _DefStatistic
variance: _DefStatistic
standard_deviation: _DefStatistic
median: _DefStatistic

#
@type_check_only
class _DefExtreme(Protocol):
    @overload
    def __call__(
        self, /, input: onp.CanArrayND[_SCT], labels: onp.ToInt | onp.ToIntND | None = None, index: None = None
    ) -> _SCT: ...
    @overload
    def __call__(
        self, /, input: onp.CanArrayND[_SCT], labels: onp.ToInt | onp.ToIntND, index: onp.ToIntND
    ) -> onp.ArrayND[_SCT]: ...
    @overload
    def __call__(
        self, /, input: onp.CanArrayND[_SCT], labels: onp.ToInt | onp.ToIntND | None = None, *, index: onp.ToIntND
    ) -> onp.ArrayND[_SCT]: ...
    @overload
    def __call__(
        self, /, input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND, index: onp.ArrayND[npc.integer]
    ) -> onp.ArrayND[Any]: ...
    @overload
    def __call__(
        self,
        /,
        input: onp.ToComplex | onp.ToComplexND,
        labels: onp.ToInt | onp.ToIntND | None = None,
        *,
        index: onp.ArrayND[npc.integer],
    ) -> onp.ArrayND[Any]: ...
    @overload
    def __call__(
        self,
        /,
        input: onp.CanArrayND[_SCT],
        labels: onp.ToInt | onp.ToIntND | None = None,
        index: onp.ToInt | onp.ToIntND | None = None,
    ) -> _SCT | onp.ArrayND[_SCT]: ...
    @overload
    def __call__(
        self,
        /,
        input: onp.ToComplex | onp.ToComplexND,
        labels: onp.ToInt | onp.ToIntND | None = None,
        index: onp.ToInt | onp.ToIntND | None = None,
    ) -> Incomplete: ...

minimum: _DefExtreme
maximum: _DefExtreme

#
@overload
def extrema(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND, labels: onp.ToInt | onp.ToIntND | None = None, index: onp.ToInt | None = None
) -> _Extrema0D[np.float64]: ...
@overload
def extrema(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    labels: onp.ToInt | onp.ToIntND | None = None,
    index: onp.ToInt | None = None,
) -> _Extrema0D[np.complex128]: ...
@overload
def extrema(
    input: _ArrayLike[_SCT], labels: onp.ToInt | onp.ToIntND | None = None, index: onp.ToInt | None = None
) -> _Extrema0D[_SCT]: ...
@overload
def extrema(
    input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND | None = None, index: onp.ToInt | None = None
) -> _Extrema0D[Any]: ...
@overload
def extrema(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND, labels: onp.ToInt | onp.ToIntND, index: onp.ToIntND
) -> _ExtremaND[np.float64]: ...
@overload
def extrema(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND, labels: onp.ToInt | onp.ToIntND | None = None, *, index: onp.ToIntND
) -> _ExtremaND[np.float64]: ...
@overload
def extrema(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND, labels: onp.ToInt | onp.ToIntND, index: onp.ToIntND
) -> _ExtremaND[np.complex128]: ...
@overload
def extrema(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND, labels: onp.ToInt | onp.ToIntND | None = None, *, index: onp.ToIntND
) -> _ExtremaND[np.complex128]: ...
@overload
def extrema(input: _ArrayLike[_SCT], labels: onp.ToInt | onp.ToIntND, index: onp.ToIntND) -> _ExtremaND[_SCT]: ...
@overload
def extrema(
    input: _ArrayLike[_SCT], labels: onp.ToInt | onp.ToIntND | None = None, *, index: onp.ToIntND
) -> _ExtremaND[_SCT]: ...
@overload
def extrema(input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND, index: onp.ToIntND) -> _ExtremaND[Any]: ...
@overload
def extrema(
    input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND | None = None, *, index: onp.ToIntND
) -> _ExtremaND[Any]: ...
@overload
def extrema(
    input: _ArrayLike[_SCT], labels: onp.ToInt | onp.ToIntND | None = None, index: onp.ToInt | onp.ToIntND | None = None
) -> _Extrema0D[_SCT] | _ExtremaND[_SCT]: ...
@overload
def extrema(
    input: onp.ToComplex | onp.ToComplexND,
    labels: onp.ToInt | onp.ToIntND | None = None,
    index: onp.ToInt | onp.ToIntND | None = None,
) -> _Extrema0D[Any] | _ExtremaND[Any]: ...

#
@type_check_only
class _DefArgExtreme(Protocol):
    @overload
    def __call__(
        self, /, input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND | None = None, index: None = None
    ) -> _Idx0D: ...
    @overload
    def __call__(
        self, /, input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND, index: onp.ArrayND[npc.integer]
    ) -> _IdxND: ...
    @overload
    def __call__(
        self,
        /,
        input: onp.ToComplex | onp.ToComplexND,
        labels: onp.ToInt | onp.ToIntND | None = None,
        *,
        index: onp.ArrayND[npc.integer],
    ) -> _IdxND: ...
    @overload
    def __call__(
        self,
        /,
        input: onp.ToComplex | onp.ToComplexND,
        labels: onp.ToInt | onp.ToIntND | None = None,
        index: onp.ToInt | onp.ToIntND | None = None,
    ) -> _Idx0D | _IdxND: ...

minimum_position: _DefArgExtreme
maximum_position: _DefArgExtreme

#
@overload
def center_of_mass(
    input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND | None = None, index: onp.ToInt | None = None
) -> _Coord0D: ...
@overload
def center_of_mass(
    input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND, index: Sequence[onp.ToInt]
) -> _Coord1D: ...
@overload
def center_of_mass(
    input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND | None = None, *, index: Sequence[onp.ToInt]
) -> _Coord1D: ...
@overload
def center_of_mass(
    input: onp.ToComplex | onp.ToComplexND, labels: onp.ToInt | onp.ToIntND, index: Sequence[Sequence[onp.ToInt | onp.ToIntND]]
) -> _CoordND: ...
@overload
def center_of_mass(
    input: onp.ToComplex | onp.ToComplexND,
    labels: onp.ToInt | onp.ToIntND | None = None,
    *,
    index: Sequence[Sequence[onp.ToInt | onp.ToIntND]],
) -> _CoordND: ...
@overload
def center_of_mass(
    input: onp.ToComplex | onp.ToComplexND,
    labels: onp.ToInt | onp.ToIntND | None = None,
    index: onp.ToInt | onp.ToIntND | None = None,
) -> _Coord0D | _Coord1D | _CoordND: ...

#
@overload
def histogram(
    input: onp.ToComplex | onp.ToComplexND,
    min: onp.ToInt,
    max: onp.ToInt,
    bins: onp.ToInt,
    labels: onp.ToInt | onp.ToIntND | None = None,
    index: onp.ToInt | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def histogram(
    input: onp.ToComplex | onp.ToComplexND,
    min: onp.ToInt,
    max: onp.ToInt,
    bins: onp.ToInt,
    labels: onp.ToInt | onp.ToIntND,
    index: onp.ToIntND,
) -> onp.ArrayND[np.object_]: ...
@overload
def histogram(
    input: onp.ToComplex | onp.ToComplexND,
    min: onp.ToInt,
    max: onp.ToInt,
    bins: onp.ToInt,
    labels: onp.ToInt | onp.ToIntND | None = None,
    *,
    index: onp.ToIntND,
) -> onp.ArrayND[np.object_]: ...
@overload
def histogram(
    input: onp.ToComplex | onp.ToComplexND,
    min: onp.ToInt,
    max: onp.ToInt,
    bins: onp.ToInt,
    labels: onp.ToInt | onp.ToIntND | None,
    index: onp.ToInt | onp.ToIntND | None,
) -> onp.ArrayND[Any]: ...

#
def watershed_ift(
    input: _ArrayLike[np.uint8 | np.uint16],
    markers: _ArrayLike[npc.signedinteger] | onp.SequenceND[int],
    structure: onp.ToInt | onp.ToIntND | None = None,
    output: onp.ArrayND[npc.signedinteger] | None = None,
) -> onp.ArrayND[npc.signedinteger]: ...
