import sys
from _typeshed import SupportsWrite as SupportsWrite
from collections.abc import Sequence
from typing import Any, Protocol, TypeVar, type_check_only
from typing_extensions import TypeAlias

import numpy as np
from numpy.typing import NDArray

from .lib import Geometry

if sys.version_info >= (3, 12):
    from collections.abc import Buffer

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)
_DType = TypeVar("_DType", bound=np.dtype[Any])
_DType_co = TypeVar("_DType_co", covariant=True, bound=np.dtype[Any])

GeoT = TypeVar("GeoT", bound=Geometry)  # noqa: Y001
OptGeoT = TypeVar("OptGeoT", bound=Geometry | None)  # noqa: Y001

@type_check_only
class SupportsArray(Protocol[_DType_co]):
    def __array__(self) -> np.ndarray[Any, _DType_co]: ...

# TODO: revisit when mypy is happy with generic recursive type alias
# NestedSequence: TypeAlias = Sequence[_T] | Sequence[NestedSequence[_T]]
NestedSequence: TypeAlias = Sequence[_T] | Sequence[Sequence[_T]] | Sequence[Sequence[Sequence[_T]]]
DualArrayLike: TypeAlias = SupportsArray[_DType] | NestedSequence[SupportsArray[_DType]] | NestedSequence[_T]

# array-like sequences: objects accepted by np.array that produce at least 1-D arrays
if sys.version_info >= (3, 12):
    ArrayLikeSeq: TypeAlias = Buffer | DualArrayLike[np.dtype[Any], _T]
else:
    ArrayLikeSeq: TypeAlias = DualArrayLike[np.dtype[Any], _T]
GeoArrayLikeSeq: TypeAlias = ArrayLikeSeq[Geometry]
OptGeoArrayLikeSeq: TypeAlias = ArrayLikeSeq[Geometry | None]

# array-like: objects accepted by np.array that may also produce 0-D array
ArrayLike: TypeAlias = _T | ArrayLikeSeq[_T]
GeoArrayLike: TypeAlias = ArrayLike[Geometry]
OptGeoArrayLike: TypeAlias = ArrayLike[Geometry | None]

# There is no way to pronounce "array of BaseGeometry" currently because of the restriction on
# NDArray type variable to np.dtype and because np.object_ is not generic.
# Note the use of `BaseGeometry` instead of `Geometry` as the alias is used in return types.
GeoArray: TypeAlias = NDArray[np.object_]

@type_check_only
class SupportsGeoInterface(Protocol):
    @property
    def __geo_interface__(self) -> dict[str, Any]: ...

# Unlike _typeshed.SupportsRead, this protocol does not require a length parameter
@type_check_only
class SupportsRead(Protocol[_T_co]):
    def read(self) -> _T_co: ...
