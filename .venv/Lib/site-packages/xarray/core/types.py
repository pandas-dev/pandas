from __future__ import annotations

import datetime
import sys
from collections.abc import Callable, Collection, Hashable, Iterator, Mapping, Sequence
from types import EllipsisType
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Protocol,
    SupportsIndex,
    TypeVar,
    Union,
)

import numpy as np
import pandas as pd

try:
    if sys.version_info >= (3, 11):
        from typing import Self, TypeAlias
    else:
        from typing import TypeAlias

        from typing_extensions import Self
except ImportError:
    if TYPE_CHECKING:
        raise
    else:
        Self: Any = None


from numpy._typing import _SupportsDType
from numpy.typing import ArrayLike

if TYPE_CHECKING:
    from xarray.backends.common import BackendEntrypoint
    from xarray.core.alignment import Aligner
    from xarray.core.common import AbstractArray, DataWithCoords
    from xarray.core.coordinates import Coordinates
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.indexes import Index, Indexes
    from xarray.core.utils import Frozen
    from xarray.core.variable import IndexVariable, Variable
    from xarray.groupers import Grouper, TimeResampler

    GroupInput: TypeAlias = (
        str
        | DataArray
        | IndexVariable
        | Sequence[Hashable]
        | Mapping[Any, Grouper]
        | None
    )

    try:
        from dask.array import Array as DaskArray
    except ImportError:
        DaskArray = np.ndarray  # type: ignore[misc, assignment, unused-ignore]

    try:
        from cubed import Array as CubedArray
    except ImportError:
        CubedArray = np.ndarray

    try:
        from zarr.core import Array as ZarrArray
    except ImportError:
        ZarrArray = np.ndarray

    # Anything that can be coerced to a shape tuple
    _ShapeLike = Union[SupportsIndex, Sequence[SupportsIndex]]
    _DTypeLikeNested = Any  # TODO: wait for support for recursive types

    # Xarray requires a Mapping[Hashable, dtype] in many places which
    # conflicts with numpys own DTypeLike (with dtypes for fields).
    # https://numpy.org/devdocs/reference/typing.html#numpy.typing.DTypeLike
    # This is a copy of this DTypeLike that allows only non-Mapping dtypes.
    DTypeLikeSave = Union[
        np.dtype[Any],
        # default data type (float64)
        None,
        # array-scalar types and generic types
        type[Any],
        # character codes, type strings or comma-separated fields, e.g., 'float64'
        str,
        # (flexible_dtype, itemsize)
        tuple[_DTypeLikeNested, int],
        # (fixed_dtype, shape)
        tuple[_DTypeLikeNested, _ShapeLike],
        # (base_dtype, new_dtype)
        tuple[_DTypeLikeNested, _DTypeLikeNested],
        # because numpy does the same?
        list[Any],
        # anything with a dtype attribute
        _SupportsDType[np.dtype[Any]],
    ]

else:
    DTypeLikeSave: Any = None

# https://mypy.readthedocs.io/en/stable/common_issues.html#variables-vs-type-aliases
try:
    from cftime import datetime as CFTimeDatetime
except ImportError:
    CFTimeDatetime = np.datetime64

DatetimeLike: TypeAlias = (
    pd.Timestamp | datetime.datetime | np.datetime64 | CFTimeDatetime
)


class Alignable(Protocol):
    """Represents any Xarray type that supports alignment.

    It may be ``Dataset``, ``DataArray`` or ``Coordinates``. This protocol class
    is needed since those types do not all have a common base class.

    """

    @property
    def dims(self) -> Frozen[Hashable, int] | tuple[Hashable, ...]: ...

    @property
    def sizes(self) -> Mapping[Hashable, int]: ...

    @property
    def xindexes(self) -> Indexes[Index]: ...

    def _reindex_callback(
        self,
        aligner: Aligner,
        dim_pos_indexers: dict[Hashable, Any],
        variables: dict[Hashable, Variable],
        indexes: dict[Hashable, Index],
        fill_value: Any,
        exclude_dims: frozenset[Hashable],
        exclude_vars: frozenset[Hashable],
    ) -> Self: ...

    def _overwrite_indexes(
        self,
        indexes: Mapping[Any, Index],
        variables: Mapping[Any, Variable] | None = None,
    ) -> Self: ...

    def __len__(self) -> int: ...

    def __iter__(self) -> Iterator[Hashable]: ...

    def copy(
        self,
        deep: bool = False,
    ) -> Self: ...


T_Alignable = TypeVar("T_Alignable", bound="Alignable")

T_Backend = TypeVar("T_Backend", bound="BackendEntrypoint")
T_Dataset = TypeVar("T_Dataset", bound="Dataset")
T_DataArray = TypeVar("T_DataArray", bound="DataArray")
T_Variable = TypeVar("T_Variable", bound="Variable")
T_Coordinates = TypeVar("T_Coordinates", bound="Coordinates")
T_Array = TypeVar("T_Array", bound="AbstractArray")
T_Index = TypeVar("T_Index", bound="Index")

# `T_Xarray` is a type variable that can be either "DataArray" or "Dataset". When used
# in a function definition, all inputs and outputs annotated with `T_Xarray` must be of
# the same concrete type, either "DataArray" or "Dataset". This is generally preferred
# over `T_DataArrayOrSet`, given the type system can determine the exact type.
T_Xarray = TypeVar("T_Xarray", "DataArray", "Dataset")

# `T_DataArrayOrSet` is a type variable that is bounded to either "DataArray" or
# "Dataset". Use it for functions that might return either type, but where the exact
# type cannot be determined statically using the type system.
T_DataArrayOrSet = TypeVar("T_DataArrayOrSet", bound=Union["Dataset", "DataArray"])

# For working directly with `DataWithCoords`. It will only allow using methods defined
# on `DataWithCoords`.
T_DataWithCoords = TypeVar("T_DataWithCoords", bound="DataWithCoords")


# Temporary placeholder for indicating an array api compliant type.
# hopefully in the future we can narrow this down more:
T_DuckArray = TypeVar("T_DuckArray", bound=Any, covariant=True)

# For typing pandas extension arrays.
T_ExtensionArray = TypeVar("T_ExtensionArray", bound=pd.api.extensions.ExtensionArray)


ScalarOrArray = Union["ArrayLike", np.generic]
VarCompatible = Union["Variable", "ScalarOrArray"]
DaCompatible = Union["DataArray", "VarCompatible"]
DsCompatible = Union["Dataset", "DaCompatible"]
GroupByCompatible = Union["Dataset", "DataArray"]

# Don't change to Hashable | Collection[Hashable]
# Read: https://github.com/pydata/xarray/issues/6142
Dims = Union[str, Collection[Hashable], EllipsisType, None]

# FYI in some cases we don't allow `None`, which this doesn't take account of.
# FYI the `str` is for a size string, e.g. "16MB", supported by dask.
T_ChunkDim: TypeAlias = str | int | Literal["auto"] | None | tuple[int, ...]
T_ChunkDimFreq: TypeAlias = Union["TimeResampler", T_ChunkDim]
T_ChunksFreq: TypeAlias = T_ChunkDim | Mapping[Any, T_ChunkDimFreq]
# We allow the tuple form of this (though arguably we could transition to named dims only)
T_Chunks: TypeAlias = T_ChunkDim | Mapping[Any, T_ChunkDim]
T_NormalizedChunks = tuple[tuple[int, ...], ...]

DataVars = Mapping[Any, Any]


ErrorOptions = Literal["raise", "ignore"]
ErrorOptionsWithWarn = Literal["raise", "warn", "ignore"]

CompatOptions = Literal[
    "identical", "equals", "broadcast_equals", "no_conflicts", "override", "minimal"
]
ConcatOptions = Literal["all", "minimal", "different"]
CombineAttrsOptions = Union[
    Literal["drop", "identical", "no_conflicts", "drop_conflicts", "override"],
    Callable[..., Any],
]
JoinOptions = Literal["outer", "inner", "left", "right", "exact", "override"]

Interp1dOptions = Literal[
    "linear", "nearest", "zero", "slinear", "quadratic", "cubic", "polynomial"
]
InterpolantOptions = Literal["barycentric", "krogh", "pchip", "spline", "akima"]
InterpOptions = Union[Interp1dOptions, InterpolantOptions]

DatetimeUnitOptions = Literal[
    "Y", "M", "W", "D", "h", "m", "s", "ms", "us", "μs", "ns", "ps", "fs", "as", None
]
NPDatetimeUnitOptions = Literal["D", "h", "m", "s", "ms", "us", "ns"]

QueryEngineOptions = Literal["python", "numexpr", None]
QueryParserOptions = Literal["pandas", "python"]

ReindexMethodOptions = Literal["nearest", "pad", "ffill", "backfill", "bfill", None]

PadModeOptions = Literal[
    "constant",
    "edge",
    "linear_ramp",
    "maximum",
    "mean",
    "median",
    "minimum",
    "reflect",
    "symmetric",
    "wrap",
]
T_PadConstantValues = float | tuple[float, float]
T_VarPadConstantValues = T_PadConstantValues | Mapping[Any, T_PadConstantValues]
T_DatasetPadConstantValues = (
    T_VarPadConstantValues | Mapping[Any, T_VarPadConstantValues]
)
PadReflectOptions = Literal["even", "odd", None]

CFCalendar = Literal[
    "standard",
    "gregorian",
    "proleptic_gregorian",
    "noleap",
    "365_day",
    "360_day",
    "julian",
    "all_leap",
    "366_day",
]

CoarsenBoundaryOptions = Literal["exact", "trim", "pad"]
SideOptions = Literal["left", "right"]
InclusiveOptions = Literal["both", "neither", "left", "right"]

ScaleOptions = Literal["linear", "symlog", "log", "logit", None]
HueStyleOptions = Literal["continuous", "discrete", None]
AspectOptions = Union[Literal["auto", "equal"], float, None]
ExtendOptions = Literal["neither", "both", "min", "max", None]

# TODO: Wait until mypy supports recursive objects in combination with typevars
_T = TypeVar("_T")
NestedSequence = Union[
    _T,
    Sequence[_T],
    Sequence[Sequence[_T]],
    Sequence[Sequence[Sequence[_T]]],
    Sequence[Sequence[Sequence[Sequence[_T]]]],
]


QuantileMethods = Literal[
    "inverted_cdf",
    "averaged_inverted_cdf",
    "closest_observation",
    "interpolated_inverted_cdf",
    "hazen",
    "weibull",
    "linear",
    "median_unbiased",
    "normal_unbiased",
    "lower",
    "higher",
    "midpoint",
    "nearest",
]


NetcdfWriteModes = Literal["w", "a"]
ZarrWriteModes = Literal["w", "w-", "a", "a-", "r+", "r"]

GroupKey = Any
GroupIndex = Union[int, slice, list[int]]
GroupIndices = tuple[GroupIndex, ...]
Bins = Union[
    int, Sequence[int], Sequence[float], Sequence[pd.Timestamp], np.ndarray, pd.Index
]

ResampleCompatible: TypeAlias = str | datetime.timedelta | pd.Timedelta | pd.DateOffset
