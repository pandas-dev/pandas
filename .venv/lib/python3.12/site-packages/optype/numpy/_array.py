import sys
from collections.abc import Mapping
from typing import Any, Protocol, Self, final

if sys.version_info >= (3, 13):
    from typing import TypeAliasType, TypeVar, runtime_checkable
else:
    from typing_extensions import TypeAliasType, TypeVar, runtime_checkable

import numpy as np

from ._shape import AnyShape, Shape
from optype._utils import set_module

__all__ = [  # noqa: RUF022
    "Array", "Array0D", "Array1D", "Array1D", "Array2D", "Array3D", "ArrayND",
    "MArray", "MArray0D", "MArray1D", "MArray1D", "MArray2D", "MArray3D",
    "Matrix",
    "CanArray", "CanArray0",
    "CanArrayND", "CanArray0D", "CanArray1D", "CanArray2D", "CanArray3D",
    "CanArrayFinalize", "CanArrayWrap",
    "HasArrayInterface", "HasArrayPriority",
]  # fmt: skip


def __dir__() -> list[str]:
    return __all__


###

_NDT = TypeVar("_NDT", bound=Shape, default=AnyShape)
_NDT_co = TypeVar("_NDT_co", bound=Shape, default=AnyShape, covariant=True)
_DTT = TypeVar("_DTT", bound=np.dtype[Any], default=np.dtype[Any])
_DTT_co = TypeVar("_DTT_co", bound=np.dtype[Any], default=np.dtype[Any], covariant=True)
_SCT = TypeVar("_SCT", bound=np.generic, default=Any)
_SCT_co = TypeVar("_SCT_co", bound=np.generic, default=Any, covariant=True)
_SCT0_co = TypeVar("_SCT0_co", bound=np.generic, covariant=True)


Matrix = TypeAliasType(
    "Matrix",
    np.matrix[tuple[int, int], np.dtype[_SCT]],
    type_params=(_SCT,),
)
"""
Alias of `np.matrix` that is similar to `ArrayND`:

```py
type Matrix[
    SCT: np.generic = Any,
] = np.matrix[tuple[int, int], np.dtype[SCT]]
```
"""

Array = TypeAliasType(
    "Array",
    np.ndarray[_NDT, np.dtype[_SCT]],
    type_params=(_NDT, _SCT),
)
"""
Shape-typed array alias, defined as:

```py
type Array[
    NDT: tuple[int, ...] = tuple[Any, ...],
    SCT: np.generic = Any,
] = np.ndarray[NDT, np.dtype[SCT]]
```
"""

ArrayND = TypeAliasType(
    "ArrayND",
    np.ndarray[_NDT, np.dtype[_SCT]],
    type_params=(_SCT, _NDT),
)
"""
Like `Array`, but with flipped type-parameters, i.e.:

type ArrayND[
    SCT: np.generic = Any,
    NDT: tuple[int, ...] = tuple[Any, ...],
] = np.ndarray[NDT, np.dtype[SCT]]

Because the optional shape-type parameter comes *after* the scalar-type, `ArrayND`
can be seen as a flexible generalization of `npt.NDArray`.
"""

MArray = TypeAliasType(
    "MArray",
    np.ma.MaskedArray[_NDT, np.dtype[_SCT]],
    type_params=(_SCT, _NDT),
)
"""
Just like `ArrayND`, but for `np.ma.MaskedArray` instead of `np.ndarray`.

type MArray[
    SCT: np.generic = Any,
    NDT: (int, ...) = (Any, ...),
] = np.ma.MaskedArray[NDT, np.dtype[SCT]]
"""


# NOTE: Before NumPy 2.1 the shape type parameter of `numpy.ndarray` was invariant. This
# lead to various issues, so we ignore that, and suppress two pyright errors that are
# reported when `numpy<2.1` is installed (but inline `# pyright: ignore` and
# `# type: ignore[misc]` will be reported as unused ignores on other versions).

# mypy: disable-error-code=misc
# pyright: reportInvalidTypeVarUse=false


@runtime_checkable
@set_module("optype.numpy")
class CanArray(Protocol[_NDT_co, _DTT_co]):
    def __array__(self, /) -> np.ndarray[_NDT_co, _DTT_co]: ...


@runtime_checkable
@set_module("optype.numpy")
class CanArrayND(Protocol[_SCT0_co, _NDT_co]):
    """
    Similar to `onp.CanArray`, but must be sized (i.e. excludes scalars), and is
    parameterized by only the scalar type (instead of the shape and dtype).
    """

    def __len__(self, /) -> int: ...
    def __array__(self, /) -> np.ndarray[_NDT_co, np.dtype[_SCT0_co]]: ...


Array0D = TypeAliasType(
    "Array0D",
    np.ndarray[tuple[()], np.dtype[_SCT]],
    type_params=(_SCT,),
)
Array1D = TypeAliasType(
    "Array1D",
    np.ndarray[tuple[int], np.dtype[_SCT]],
    type_params=(_SCT,),
)
Array2D = TypeAliasType(
    "Array2D",
    np.ndarray[tuple[int, int], np.dtype[_SCT]],
    type_params=(_SCT,),
)
Array3D = TypeAliasType(
    "Array3D",
    np.ndarray[tuple[int, int, int], np.dtype[_SCT]],
    type_params=(_SCT,),
)

MArray0D = TypeAliasType(
    "MArray0D",
    np.ma.MaskedArray[tuple[()], np.dtype[_SCT]],
    type_params=(_SCT,),
)
MArray1D = TypeAliasType(
    "MArray1D",
    np.ma.MaskedArray[tuple[int], np.dtype[_SCT]],
    type_params=(_SCT,),
)
MArray2D = TypeAliasType(
    "MArray2D",
    np.ma.MaskedArray[tuple[int, int], np.dtype[_SCT]],
    type_params=(_SCT,),
)
MArray3D = TypeAliasType(
    "MArray3D",
    np.ma.MaskedArray[tuple[int, int, int], np.dtype[_SCT]],
    type_params=(_SCT,),
)


###########################


@final
@runtime_checkable
@set_module("optype.numpy")
class CanArray0(Protocol[_SCT_co]):
    """
    The 0-d variant of `optype.numpy.CanArray`.

    This accepts e.g. `np.array(4)` and `np.int8(2)`, but rejects `np.array([4, 2])`.

    Note that this is different from `optype.numpy.CanArray0D`, which would reject
    "bare" scalars like `np.int8(2)`.
    """

    def __array__(self, /) -> np.ndarray[tuple[()], np.dtype[_SCT_co]]: ...


@runtime_checkable
@set_module("optype.numpy")
class CanArray0D(Protocol[_SCT_co]):
    """
    The 0-d variant of `optype.numpy.CanArrayND`.

    This accepts e.g. `np.array(3.14)`, but rejects `np.float64(3.14)` .

    Note that this is different from `optype.numpy.CanArray0`, which would accept
    "bare" scalars like `np.float64(3.14)`.
    """

    def __len__(self, /) -> int: ...  # always 0
    def __array__(self, /) -> np.ndarray[tuple[()], np.dtype[_SCT_co]]: ...


@runtime_checkable
@set_module("optype.numpy")
class CanArray1D(Protocol[_SCT_co]):
    """The 1-d variant of `optype.numpy.CanArrayND`."""

    def __len__(self, /) -> int: ...
    def __array__(self, /) -> np.ndarray[tuple[int], np.dtype[_SCT_co]]: ...


@runtime_checkable
@set_module("optype.numpy")
class CanArray2D(Protocol[_SCT_co]):
    """The 2-d variant of `optype.numpy.CanArrayND`."""

    def __len__(self, /) -> int: ...
    def __array__(self, /) -> np.ndarray[tuple[int, int], np.dtype[_SCT_co]]: ...


@runtime_checkable
@set_module("optype.numpy")
class CanArray3D(Protocol[_SCT_co]):
    """The 2-d variant of `optype.numpy.CanArrayND`."""

    def __len__(self, /) -> int: ...
    def __array__(self, /) -> np.ndarray[tuple[int, int, int], np.dtype[_SCT_co]]: ...


# this is almost always a `ndarray`, but setting a `bound` might break in some
# edge cases
_T_contra = TypeVar("_T_contra", contravariant=True, default=Any)


@runtime_checkable
@set_module("optype.numpy")
class CanArrayFinalize(Protocol[_T_contra]):
    def __array_finalize__(self, obj: _T_contra, /) -> None: ...


@runtime_checkable
@set_module("optype.numpy")
class CanArrayWrap(Protocol):
    def __array_wrap__(
        self,
        array: np.ndarray[_NDT, _DTT],
        context: tuple[np.ufunc, tuple[Any, ...], int] | None = ...,
        /,
    ) -> np.ndarray[_NDT, _DTT] | Self: ...


_ArrayInterfaceT_co = TypeVar(
    "_ArrayInterfaceT_co",
    covariant=True,
    bound=Mapping[str, object],
    default=dict[str, Any],
)


@runtime_checkable
@set_module("optype.numpy")
class HasArrayInterface(Protocol[_ArrayInterfaceT_co]):
    @property
    def __array_interface__(self, /) -> _ArrayInterfaceT_co: ...


@runtime_checkable
@set_module("optype.numpy")
class HasArrayPriority(Protocol):
    @property
    def __array_priority__(self, /) -> float: ...
