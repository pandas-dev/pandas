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
from types import EllipsisType
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    ParamSpec,
    Protocol,
    SupportsIndex,
    TypeVar,
    overload,
)

import numpy as np
import numpy.typing as npt

# To prevent import cycles place any internal imports in the branch below
# and use a string literal forward reference to it in subsequent types
# https://mypy.readthedocs.io/en/stable/runtime_troubles.html#import-cycles

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

    type ScalarLike_co = int | float | complex | str | bytes | np.generic

    # numpy compatible types
    type NumpyValueArrayLike = ScalarLike_co | npt.ArrayLike
    type NumpySorter = npt.NDArray[np.integer] | None


P = ParamSpec("P")

HashableT = TypeVar("HashableT", bound=Hashable)
HashableT2 = TypeVar("HashableT2", bound=Hashable)
MutableMappingT = TypeVar("MutableMappingT", bound=MutableMapping[Any, Any])

# array-like

type ArrayLike = ExtensionArray | np.ndarray
ArrayLikeT = TypeVar("ArrayLikeT", "ExtensionArray", np.ndarray)
type AnyArrayLike = ArrayLike | Index | Series
type TimeArrayLike = DatetimeArray | TimedeltaArray

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


type ListLike = AnyArrayLike | SequenceNotStr[object] | range

# scalars

type PythonScalar = str | float | bool
type DatetimeLikeScalar = Period | Timestamp | Timedelta

# aligned with pandas-stubs - typical scalars found in Series.  Explicitly leaves
# out object
type _IndexIterScalar = (
    str
    | bytes
    | date
    | datetime
    | timedelta
    | np.datetime64
    | np.timedelta64
    | bool
    | int
    | float
    | Timestamp
    | Timedelta
)
type Scalar = (
    _IndexIterScalar
    | Interval[Any]
    | complex
    | np.integer
    | np.floating
    | np.complexfloating
)

IntStrT = TypeVar("IntStrT", bound=int | str)

# timestamp and timedelta convertible types

type TimestampConvertibleTypes = (
    Timestamp | date | np.datetime64 | np.int64 | float | str
)
type TimestampNonexistent = (
    Literal["shift_forward", "shift_backward", "NaT", "raise"] | timedelta
)

type TimedeltaConvertibleTypes = (
    Timedelta | timedelta | np.timedelta64 | np.int64 | float | str
)
type Timezone = str | tzinfo

type ToTimestampHow = Literal["s", "e", "start", "end"]

# NDFrameT is stricter and ensures that the same subclass of NDFrame always is
# used. E.g. `def func(a: NDFrameT) -> NDFrameT: ...` means that if a
# Series is passed into a function, a Series is always returned and if a DataFrame is
# passed in, a DataFrame is always returned.
NDFrameT = TypeVar("NDFrameT", bound="NDFrame")

IndexT = TypeVar("IndexT", bound="Index")
FreqIndexT = TypeVar("FreqIndexT", "DatetimeIndex", "PeriodIndex", "TimedeltaIndex")
NumpyIndexT = TypeVar("NumpyIndexT", np.ndarray, "Index")

type AxisInt = int
type Axis = AxisInt | Literal["index", "columns", "rows"]
type IndexLabel = Hashable | Sequence[Hashable]
type Level = Hashable
type Shape = tuple[int, ...]
type Suffixes = SequenceNotStr[str | None]
type Ordered = bool | None
type JSONSerializable = PythonScalar | list[Any] | dict[Any, Any] | None
type Frequency = str | BaseOffset
type Axes = ListLike

type RandomState = (
    int
    | np.ndarray
    | np.random.Generator
    | np.random.BitGenerator
    | np.random.RandomState
)


# dtypes
type NpDtype = str | np.dtype | type[str | complex | bool | object]
type Dtype = ExtensionDtype | NpDtype
type AstypeArg = ExtensionDtype | npt.DTypeLike
# DtypeArg specifies all allowable dtypes in a functions its dtype argument
type DtypeArg = Dtype | Mapping[Hashable, Dtype]
type DtypeObj = np.dtype | ExtensionDtype

# converters
type ConvertersArg = dict[Hashable, Callable[[Dtype], Dtype]]

# parse_dates
type ParseDatesArg = (
    bool | list[Hashable] | list[list[Hashable]] | dict[Hashable, list[Hashable]]
)

# For functions like rename that convert one label to another
type Renamer = Mapping[Any, Hashable] | Callable[[Any], Hashable]

# to maintain type information across generic functions and parametrization
T = TypeVar("T")

# used in decorators to preserve the signature of the function it decorates
# see https://mypy.readthedocs.io/en/stable/generics.html#declaring-decorators
type FuncType = Callable[..., Any]
F = TypeVar("F", bound=FuncType)
TypeT = TypeVar("TypeT", bound=type)

# types of vectorized key functions for DataFrame::sort_values and
# DataFrame::sort_index, among others
type ValueKeyFunc = Callable[[Series], Series | AnyArrayLike] | None
type IndexKeyFunc = Callable[[Index], Index | AnyArrayLike] | None

# types of `func` kwarg for DataFrame.aggregate and Series.aggregate
type AggFuncTypeBase = Callable[..., Any] | str
type AggFuncTypeDict = MutableMapping[Hashable, AggFuncTypeBase | list[AggFuncTypeBase]]
type AggFuncType = AggFuncTypeBase | list[AggFuncTypeBase] | AggFuncTypeDict
type AggObjType = (
    Series
    | DataFrame
    | GroupBy[Any]
    | SeriesGroupBy
    | DataFrameGroupBy
    | BaseWindow
    | Resampler
)

type PythonFuncType = Callable[[Any], Any]

# filenames and file-like-objects
AnyStr_co = TypeVar("AnyStr_co", str, bytes, covariant=True)
AnyStr_contra = TypeVar("AnyStr_contra", str, bytes, contravariant=True)


class BaseBuffer(Protocol):
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


type FilePath = str | PathLike[str]

# for arbitrary kwargs passed during reading/writing files
type StorageOptions = dict[str, Any] | None

# compression keywords and compression
type CompressionDict = dict[str, Any]
type CompressionOptions = (
    Literal["infer", "gzip", "bz2", "zip", "xz", "zstd", "tar"] | CompressionDict | None
)
type ParquetCompressionOptions = (
    Literal["snappy", "gzip", "brotli", "lz4", "zstd"] | None
)

# types in DataFrameFormatter
type FormattersType = (
    list[Callable[..., Any]]
    | tuple[Callable[..., Any], ...]
    | Mapping[str | int, Callable[..., Any]]
)
type ColspaceType = Mapping[Hashable, str | int]
type FloatFormatType = str | Callable[..., Any] | EngFormatter
type ColspaceArgType = str | int | Sequence[str | int] | Mapping[Hashable, str | int]

# Arguments for fillna()
type FillnaOptions = Literal["backfill", "bfill", "ffill", "pad"]
type InterpolateOptions = Literal[
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
type Manager = BlockManager | SingleBlockManager

# indexing
# PositionalIndexer -> valid 1D positional indexer, e.g. can pass
# to ndarray.__getitem__
# ScalarIndexer is for a single value as the index
# SequenceIndexer is for list like or slices (but not tuples)
# PositionalIndexerTuple is extends the PositionalIndexer for 2D arrays
# These are used in various __getitem__ overloads
# Using List[int] here rather than Sequence[int] to disallow tuples.
type ScalarIndexer = int | np.integer
type SequenceIndexer = slice | EllipsisType | list[int] | np.ndarray
type PositionalIndexer = ScalarIndexer | SequenceIndexer
type PositionalIndexerTuple = tuple[PositionalIndexer, PositionalIndexer]
type PositionalIndexer2D = PositionalIndexer | PositionalIndexerTuple
type TakeIndexer = Sequence[int] | Sequence[np.integer] | npt.NDArray[np.integer]

# Shared by functions such as drop and astype
type IgnoreRaise = Literal["ignore", "raise"]

# Windowing rank methods
type WindowingRankType = Literal["average", "min", "max"]

# rank tie-breaking methods (Series.rank / DataFrame.rank / GroupBy.rank)
type RankMethod = Literal["average", "min", "max", "first", "dense"]
type RankNaOption = Literal["keep", "top", "bottom"]

# read_csv engines
type CSVEngine = Literal["c", "python", "pyarrow", "python-fwf"]

# read_json engines
type JSONEngine = Literal["ujson", "pyarrow"]

# read_xml parsers
type XMLParsers = Literal["lxml", "etree"]

# read_html flavors
type HTMLFlavors = Literal["lxml", "html5lib", "bs4"]

# Interval closed type
type IntervalLeftRight = Literal["left", "right"]
type IntervalClosedType = IntervalLeftRight | Literal["both", "neither"]

# datetime and NaTType
type DatetimeNaTType = datetime | NaTType
type DateTimeErrorChoices = Literal["raise", "coerce"]

# sort_index
type SortKind = Literal["quicksort", "mergesort", "heapsort", "stable"]
type NaPosition = Literal["first", "last"]

# Arguments for nsmallest and nlargest
type NsmallestNlargestKeep = Literal["first", "last", "all"]

# quantile interpolation
type QuantileInterpolation = Literal["linear", "lower", "higher", "midpoint", "nearest"]

# plotting
type PlottingOrientation = Literal["horizontal", "vertical"]

# dropna
type AnyAll = Literal["any", "all"]

# merge
type MergeHow = Literal[
    "left", "right", "inner", "outer", "cross", "left_anti", "right_anti"
]
type MergeValidate = Literal[
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
type JoinHow = Literal["left", "right", "inner", "outer"]
type JoinValidate = Literal[
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
type ReindexMethod = FillnaOptions | Literal["nearest"]

type MatplotlibColor = str | Sequence[float]
type TimeGrouperOrigin = (
    Timestamp | Literal["epoch", "start", "start_day", "end", "end_day"]
)
type TimeAmbiguous = Literal["infer", "NaT", "raise"] | bool | npt.NDArray[np.bool_]
type TimeNonexistent = (
    Literal["shift_forward", "shift_backward", "NaT", "raise"] | timedelta
)

type DropKeep = Literal["first", "last", False]
type CorrelationMethod = (
    Literal["pearson", "kendall", "spearman"]
    | Callable[[np.ndarray, np.ndarray], float]
)

type AlignJoin = Literal["outer", "inner", "left", "right"]
type DtypeBackend = Literal["pyarrow", "numpy_nullable"]

type TimeUnit = Literal["s", "ms", "us", "ns"]

type NaTStrings = Literal["", "NaT", "nat", "NAT", "nan", "NaN", "NAN"]
type OpenFileErrors = Literal[
    "strict",
    "ignore",
    "replace",
    "surrogateescape",
    "xmlcharrefreplace",
    "backslashreplace",
    "namereplace",
]

# update
type UpdateJoin = Literal["left"]

# applymap
type NaAction = Literal["ignore"]

# from_dict
type FromDictOrient = Literal["columns", "index", "tight"]

# to_stata
type ToStataByteorder = Literal[">", "<", "little", "big"]

# ExcelWriter
type ExcelWriterIfSheetExists = Literal["error", "new", "replace", "overlay"]
type ExcelWriterMergeCells = bool | Literal["columns"]

# Offsets
type OffsetCalendar = np.busdaycalendar | AbstractHolidayCalendar

# read_csv: usecols
type UsecolsArgType[HashableT: Hashable] = (
    SequenceNotStr[Hashable] | range | AnyArrayLike | Callable[[HashableT], bool] | None
)

# maintain the sub-type of any hashable sequence
SequenceT = TypeVar("SequenceT", bound=Sequence[Hashable])

type SliceType = Hashable | None


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
