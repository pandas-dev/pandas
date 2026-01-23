from __future__ import annotations

from builtins import type as type_t  # pyright: ignore[reportUnusedImport]
from collections.abc import (
    Callable,
    Hashable,
    Iterator,
    Mapping,
    MutableMapping,
    Sequence,
)
from datetime import (
    date,
    datetime,
    timedelta,
    tzinfo,
)
from os import PathLike
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    ParamSpec,
    Protocol,
    SupportsIndex,
    TypeAlias,
    TypeVar,
    Union,
    overload,
)

import numpy as np
import numpy.typing as npt

# To prevent import cycles place any internal imports in the branch below
# and use a string literal forward reference to it in subsequent types
# https://mypy.readthedocs.io/en/latest/common_issues.html#import-cycles

# Note that Union is needed when a Union includes a pandas type

if TYPE_CHECKING:
    from pandas._libs import (
        NaTType,
        Period,
        Timedelta,
        Timestamp,
    )
    from pandas._libs.tslibs import BaseOffset

    from pandas.core.dtypes.dtypes import ExtensionDtype

    from pandas import (
        DatetimeIndex,
        Interval,
        PeriodIndex,
        TimedeltaIndex,
    )
    from pandas.arrays import (
        DatetimeArray,
        TimedeltaArray,
    )
    from pandas.core.arrays.base import ExtensionArray
    from pandas.core.frame import DataFrame
    from pandas.core.generic import NDFrame
    from pandas.core.groupby.generic import (
        DataFrameGroupBy,
        GroupBy,
        SeriesGroupBy,
    )
    from pandas.core.indexes.base import Index
    from pandas.core.internals import (
        BlockManager,
        SingleBlockManager,
    )
    from pandas.core.resample import Resampler
    from pandas.core.series import Series
    from pandas.core.window.rolling import BaseWindow

    from pandas.io.formats.format import EngFormatter
    from pandas.tseries.holiday import AbstractHolidayCalendar

    ScalarLike_co: TypeAlias = int | float | complex | str | bytes | np.generic

    # numpy compatible types
    NumpyValueArrayLike: TypeAlias = ScalarLike_co | npt.ArrayLike
    NumpySorter: TypeAlias = npt._ArrayLikeInt_co | None


P = ParamSpec("P")

HashableT = TypeVar("HashableT", bound=Hashable)
HashableT2 = TypeVar("HashableT2", bound=Hashable)
MutableMappingT = TypeVar("MutableMappingT", bound=MutableMapping)

# array-like

ArrayLike: TypeAlias = Union["ExtensionArray", np.ndarray]
ArrayLikeT = TypeVar("ArrayLikeT", "ExtensionArray", np.ndarray)
AnyArrayLike: TypeAlias = Union[ArrayLike, "Index", "Series"]
TimeArrayLike: TypeAlias = Union["DatetimeArray", "TimedeltaArray"]

# list-like

# from https://github.com/hauntsaninja/useful_types
# includes Sequence-like objects but excludes str and bytes
_T_co = TypeVar("_T_co", covariant=True)


class SequenceNotStr(Protocol[_T_co]):
    __module__: str = "pandas.api.typing.aliases"

    @overload
    def __getitem__(self, index: SupportsIndex, /) -> _T_co: ...

    @overload
    def __getitem__(self, index: slice, /) -> Sequence[_T_co]: ...

    def __contains__(self, value: object, /) -> bool: ...

    def __len__(self) -> int: ...

    def __iter__(self) -> Iterator[_T_co]: ...

    def index(self, value: Any, start: int = ..., stop: int = ..., /) -> int: ...

    def count(self, value: Any, /) -> int: ...

    def __reversed__(self) -> Iterator[_T_co]: ...


ListLike: TypeAlias = AnyArrayLike | SequenceNotStr | range

# scalars

PythonScalar: TypeAlias = str | float | bool
DatetimeLikeScalar: TypeAlias = Union["Period", "Timestamp", "Timedelta"]

# aligned with pandas-stubs - typical scalars found in Series.  Explicitly leaves
# out object
_IndexIterScalar: TypeAlias = Union[
    str,
    bytes,
    date,
    datetime,
    timedelta,
    np.datetime64,
    np.timedelta64,
    bool,
    int,
    float,
    "Timestamp",
    "Timedelta",
]
Scalar: TypeAlias = Union[
    _IndexIterScalar, "Interval", complex, np.integer, np.floating, np.complexfloating
]

IntStrT = TypeVar("IntStrT", bound=int | str)

# timestamp and timedelta convertible types

TimestampConvertibleTypes: TypeAlias = Union[
    "Timestamp", date, np.datetime64, np.int64, float, str
]
TimestampNonexistent: TypeAlias = (
    Literal["shift_forward", "shift_backward", "NaT", "raise"] | timedelta
)

TimedeltaConvertibleTypes: TypeAlias = Union[
    "Timedelta", timedelta, np.timedelta64, np.int64, float, str
]
Timezone: TypeAlias = str | tzinfo

ToTimestampHow: TypeAlias = Literal["s", "e", "start", "end"]

# NDFrameT is stricter and ensures that the same subclass of NDFrame always is
# used. E.g. `def func(a: NDFrameT) -> NDFrameT: ...` means that if a
# Series is passed into a function, a Series is always returned and if a DataFrame is
# passed in, a DataFrame is always returned.
NDFrameT = TypeVar("NDFrameT", bound="NDFrame")

IndexT = TypeVar("IndexT", bound="Index")
FreqIndexT = TypeVar("FreqIndexT", "DatetimeIndex", "PeriodIndex", "TimedeltaIndex")
NumpyIndexT = TypeVar("NumpyIndexT", np.ndarray, "Index")

AxisInt: TypeAlias = int
Axis: TypeAlias = AxisInt | Literal["index", "columns", "rows"]
IndexLabel: TypeAlias = Hashable | Sequence[Hashable]
Level: TypeAlias = Hashable
Shape: TypeAlias = tuple[int, ...]
Suffixes: TypeAlias = Sequence[str | None]
Ordered: TypeAlias = bool | None
JSONSerializable: TypeAlias = PythonScalar | list | dict | None
Frequency: TypeAlias = Union[str, "BaseOffset"]
Axes: TypeAlias = ListLike

RandomState: TypeAlias = (
    int
    | np.ndarray
    | np.random.Generator
    | np.random.BitGenerator
    | np.random.RandomState
)


# dtypes
NpDtype: TypeAlias = str | np.dtype | type[str | complex | bool | object]
Dtype: TypeAlias = Union["ExtensionDtype", NpDtype]
AstypeArg: TypeAlias = Union["ExtensionDtype", npt.DTypeLike]
# DtypeArg specifies all allowable dtypes in a functions its dtype argument
DtypeArg: TypeAlias = Dtype | Mapping[Hashable, Dtype]
DtypeObj: TypeAlias = Union[np.dtype, "ExtensionDtype"]

# converters
ConvertersArg: TypeAlias = dict[Hashable, Callable[[Dtype], Dtype]]

# parse_dates
ParseDatesArg: TypeAlias = (
    bool | list[Hashable] | list[list[Hashable]] | dict[Hashable, list[Hashable]]
)

# For functions like rename that convert one label to another
Renamer: TypeAlias = Mapping[Any, Hashable] | Callable[[Any], Hashable]

# to maintain type information across generic functions and parametrization
T = TypeVar("T")

# used in decorators to preserve the signature of the function it decorates
# see https://mypy.readthedocs.io/en/stable/generics.html#declaring-decorators
FuncType: TypeAlias = Callable[..., Any]
F = TypeVar("F", bound=FuncType)
TypeT = TypeVar("TypeT", bound=type)

# types of vectorized key functions for DataFrame::sort_values and
# DataFrame::sort_index, among others
ValueKeyFunc: TypeAlias = Callable[["Series"], Union["Series", AnyArrayLike]] | None
IndexKeyFunc: TypeAlias = Callable[["Index"], Union["Index", AnyArrayLike]] | None

# types of `func` kwarg for DataFrame.aggregate and Series.aggregate
AggFuncTypeBase: TypeAlias = Callable | str
AggFuncTypeDict: TypeAlias = MutableMapping[
    Hashable, AggFuncTypeBase | list[AggFuncTypeBase]
]
AggFuncType: TypeAlias = AggFuncTypeBase | list[AggFuncTypeBase] | AggFuncTypeDict
AggObjType: TypeAlias = Union[
    "Series",
    "DataFrame",
    "GroupBy",
    "SeriesGroupBy",
    "DataFrameGroupBy",
    "BaseWindow",
    "Resampler",
]

PythonFuncType: TypeAlias = Callable[[Any], Any]

# filenames and file-like-objects
AnyStr_co = TypeVar("AnyStr_co", str, bytes, covariant=True)
AnyStr_contra = TypeVar("AnyStr_contra", str, bytes, contravariant=True)


class BaseBuffer(Protocol):
    @property
    def mode(self) -> str:
        # for _get_filepath_or_buffer
        ...

    def seek(self, offset: int, whence: int = ..., /) -> int:
        # with one argument: gzip.GzipFile, bz2.BZ2File
        # with two arguments: zip.ZipFile, read_sas
        ...

    def seekable(self) -> bool:
        # for bz2.BZ2File
        ...

    def tell(self) -> int:
        # for zip.ZipFile, read_stata, to_stata
        ...


class ReadBuffer(BaseBuffer, Protocol[AnyStr_co]):
    __module__: str = "pandas.api.typing.aliases"

    def read(self, n: int = ..., /) -> AnyStr_co:
        # for BytesIOWrapper, gzip.GzipFile, bz2.BZ2File
        ...


class WriteBuffer(BaseBuffer, Protocol[AnyStr_contra]):
    __module__: str = "pandas.api.typing.aliases"

    def write(self, b: AnyStr_contra, /) -> Any:
        # for gzip.GzipFile, bz2.BZ2File
        ...

    def flush(self) -> Any:
        # for gzip.GzipFile, bz2.BZ2File
        ...


class ReadPickleBuffer(ReadBuffer[bytes], Protocol):
    __module__: str = "pandas.api.typing.aliases"

    def readline(self) -> bytes: ...


class WriteExcelBuffer(WriteBuffer[bytes], Protocol):
    __module__: str = "pandas.api.typing.aliases"

    def truncate(self, size: int | None = ..., /) -> int: ...


class ReadCsvBuffer(ReadBuffer[AnyStr_co], Protocol):
    __module__: str = "pandas.api.typing.aliases"

    def __iter__(self) -> Iterator[AnyStr_co]:
        # for engine=python
        ...

    def fileno(self) -> int:
        # for _MMapWrapper
        ...

    def readline(self) -> AnyStr_co:
        # for engine=python
        ...

    @property
    def closed(self) -> bool:
        # for engine=pyarrow
        ...


FilePath: TypeAlias = str | PathLike[str]

# for arbitrary kwargs passed during reading/writing files
StorageOptions: TypeAlias = dict[str, Any] | None

# compression keywords and compression
CompressionDict: TypeAlias = dict[str, Any]
CompressionOptions: TypeAlias = (
    Literal["infer", "gzip", "bz2", "zip", "xz", "zstd", "tar"] | CompressionDict | None
)
ParquetCompressionOptions: TypeAlias = (
    Literal["snappy", "gzip", "brotli", "lz4", "zstd"] | None
)

# types in DataFrameFormatter
FormattersType: TypeAlias = (
    list[Callable] | tuple[Callable, ...] | Mapping[str | int, Callable]
)
ColspaceType: TypeAlias = Mapping[Hashable, str | int]
FloatFormatType: TypeAlias = Union[str, Callable, "EngFormatter"]
ColspaceArgType: TypeAlias = (
    str | int | Sequence[str | int] | Mapping[Hashable, str | int]
)

# Arguments for fillna()
FillnaOptions: TypeAlias = Literal["backfill", "bfill", "ffill", "pad"]
InterpolateOptions: TypeAlias = Literal[
    "linear",
    "time",
    "index",
    "values",
    "nearest",
    "zero",
    "slinear",
    "quadratic",
    "cubic",
    "barycentric",
    "polynomial",
    "krogh",
    "piecewise_polynomial",
    "spline",
    "pchip",
    "akima",
    "cubicspline",
    "from_derivatives",
]

# internals
Manager: TypeAlias = Union["BlockManager", "SingleBlockManager"]

# indexing
# PositionalIndexer -> valid 1D positional indexer, e.g. can pass
# to ndarray.__getitem__
# ScalarIndexer is for a single value as the index
# SequenceIndexer is for list like or slices (but not tuples)
# PositionalIndexerTuple is extends the PositionalIndexer for 2D arrays
# These are used in various __getitem__ overloads
# TODO(typing#684): add Ellipsis, see
# https://github.com/python/typing/issues/684#issuecomment-548203158
# https://bugs.python.org/issue41810
# Using List[int] here rather than Sequence[int] to disallow tuples.
ScalarIndexer: TypeAlias = int | np.integer
SequenceIndexer: TypeAlias = slice | list[int] | np.ndarray
PositionalIndexer: TypeAlias = ScalarIndexer | SequenceIndexer
PositionalIndexerTuple: TypeAlias = tuple[PositionalIndexer, PositionalIndexer]
PositionalIndexer2D: TypeAlias = PositionalIndexer | PositionalIndexerTuple
TakeIndexer: TypeAlias = Sequence[int] | Sequence[np.integer] | npt.NDArray[np.integer]

# Shared by functions such as drop and astype
IgnoreRaise: TypeAlias = Literal["ignore", "raise"]

# Windowing rank methods
WindowingRankType: TypeAlias = Literal["average", "min", "max"]

# read_csv engines
CSVEngine: TypeAlias = Literal["c", "python", "pyarrow", "python-fwf"]

# read_json engines
JSONEngine: TypeAlias = Literal["ujson", "pyarrow"]

# read_xml parsers
XMLParsers: TypeAlias = Literal["lxml", "etree"]

# read_html flavors
HTMLFlavors: TypeAlias = Literal["lxml", "html5lib", "bs4"]

# Interval closed type
IntervalLeftRight: TypeAlias = Literal["left", "right"]
IntervalClosedType: TypeAlias = IntervalLeftRight | Literal["both", "neither"]

# datetime and NaTType
DatetimeNaTType: TypeAlias = Union[datetime, "NaTType"]
DateTimeErrorChoices: TypeAlias = Literal["raise", "coerce"]

# sort_index
SortKind: TypeAlias = Literal["quicksort", "mergesort", "heapsort", "stable"]
NaPosition: TypeAlias = Literal["first", "last"]

# Arguments for nsmallest and nlargest
NsmallestNlargestKeep: TypeAlias = Literal["first", "last", "all"]

# quantile interpolation
QuantileInterpolation: TypeAlias = Literal[
    "linear", "lower", "higher", "midpoint", "nearest"
]

# plotting
PlottingOrientation: TypeAlias = Literal["horizontal", "vertical"]

# dropna
AnyAll: TypeAlias = Literal["any", "all"]

# merge
MergeHow: TypeAlias = Literal[
    "left", "right", "inner", "outer", "cross", "left_anti", "right_anti"
]
MergeValidate: TypeAlias = Literal[
    "one_to_one",
    "1:1",
    "one_to_many",
    "1:m",
    "many_to_one",
    "m:1",
    "many_to_many",
    "m:m",
]

# join
JoinHow: TypeAlias = Literal["left", "right", "inner", "outer"]
JoinValidate: TypeAlias = Literal[
    "one_to_one",
    "1:1",
    "one_to_many",
    "1:m",
    "many_to_one",
    "m:1",
    "many_to_many",
    "m:m",
]

# reindex
ReindexMethod: TypeAlias = FillnaOptions | Literal["nearest"]

MatplotlibColor: TypeAlias = str | Sequence[float]
TimeGrouperOrigin: TypeAlias = Union[
    "Timestamp", Literal["epoch", "start", "start_day", "end", "end_day"]
]
TimeAmbiguous: TypeAlias = (
    Literal["infer", "NaT", "raise"] | bool | npt.NDArray[np.bool_]
)
TimeNonexistent: TypeAlias = (
    Literal["shift_forward", "shift_backward", "NaT", "raise"] | timedelta
)

DropKeep: TypeAlias = Literal["first", "last", False]
CorrelationMethod: TypeAlias = (
    Literal["pearson", "kendall", "spearman"]
    | Callable[[np.ndarray, np.ndarray], float]
)

AlignJoin: TypeAlias = Literal["outer", "inner", "left", "right"]
DtypeBackend: TypeAlias = Literal["pyarrow", "numpy_nullable"]

TimeUnit: TypeAlias = Literal["s", "ms", "us", "ns"]
OpenFileErrors: TypeAlias = Literal[
    "strict",
    "ignore",
    "replace",
    "surrogateescape",
    "xmlcharrefreplace",
    "backslashreplace",
    "namereplace",
]

# update
UpdateJoin: TypeAlias = Literal["left"]

# applymap
NaAction: TypeAlias = Literal["ignore"]

# from_dict
FromDictOrient: TypeAlias = Literal["columns", "index", "tight"]

# to_stata
ToStataByteorder: TypeAlias = Literal[">", "<", "little", "big"]

# ExcelWriter
ExcelWriterIfSheetExists: TypeAlias = Literal["error", "new", "replace", "overlay"]
ExcelWriterMergeCells: TypeAlias = bool | Literal["columns"]

# Offsets
OffsetCalendar: TypeAlias = Union[np.busdaycalendar, "AbstractHolidayCalendar"]

# read_csv: usecols
UsecolsArgType: TypeAlias = (
    SequenceNotStr[Hashable] | range | AnyArrayLike | Callable[[HashableT], bool] | None
)

# maintain the sub-type of any hashable sequence
SequenceT = TypeVar("SequenceT", bound=Sequence[Hashable])

SliceType: TypeAlias = Hashable | None


# Arrow PyCapsule Interface
# from https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html#protocol-typehints


class ArrowArrayExportable(Protocol):
    """
    An object with an ``__arrow_c_array__`` method.

    This method indicates the object is an Arrow-compatible object implementing
    the `Arrow PyCapsule Protocol`_ (exposing the `Arrow C Data Interface`_ in
    Python), enabling zero-copy Arrow data interchange across libraries.

    .. _Arrow PyCapsule Protocol: https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html
    .. _Arrow C Data Interface: https://arrow.apache.org/docs/format/CDataInterface.html

    """

    def __arrow_c_array__(
        self, requested_schema: object | None = None
    ) -> tuple[object, object]: ...


class ArrowStreamExportable(Protocol):
    """
    An object with an ``__arrow_c_stream__`` method.

    This method indicates the object is an Arrow-compatible object implementing
    the `Arrow PyCapsule Protocol`_ (exposing the `Arrow C Data Interface`_
    for streams in Python), enabling zero-copy Arrow data interchange across
    libraries.

    .. _Arrow PyCapsule Protocol: https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html
    .. _Arrow C Stream Interface: https://arrow.apache.org/docs/format/CStreamInterface.html

    """

    def __arrow_c_stream__(self, requested_schema: object | None = None) -> object: ...


__all__ = ["type_t"]
