from __future__ import annotations

import asyncio
import copy
import datetime
import io
import math
import sys
import warnings
from collections import defaultdict
from collections.abc import (
    Callable,
    Collection,
    Hashable,
    Iterable,
    Iterator,
    Mapping,
    MutableMapping,
    Sequence,
)
from functools import partial
from html import escape
from numbers import Number
from operator import methodcaller
from os import PathLike
from types import EllipsisType
from typing import IO, TYPE_CHECKING, Any, Literal, cast, overload

import numpy as np
import pandas as pd

from xarray.coding.calendar_ops import convert_calendar, interp_calendar
from xarray.coding.cftimeindex import CFTimeIndex, _parse_array_of_cftime_strings
from xarray.compat.array_api_compat import to_like_array
from xarray.computation import ops
from xarray.computation.arithmetic import DatasetArithmetic
from xarray.core import dtypes as xrdtypes
from xarray.core import duck_array_ops, formatting, formatting_html, utils
from xarray.core._aggregations import DatasetAggregations
from xarray.core.common import (
    DataWithCoords,
    _contains_datetime_like_objects,
    get_chunksizes,
)
from xarray.core.coordinates import (
    Coordinates,
    DatasetCoordinates,
    assert_coordinate_consistent,
)
from xarray.core.dataset_utils import _get_virtual_variable, _LocIndexer
from xarray.core.dataset_variables import DataVariables
from xarray.core.duck_array_ops import datetime_to_numeric
from xarray.core.indexes import (
    Index,
    Indexes,
    PandasIndex,
    PandasMultiIndex,
    assert_no_index_corrupted,
    create_default_index_implicit,
    filter_indexes_from_coords,
    isel_indexes,
    remove_unused_levels_categories,
    roll_indexes,
)
from xarray.core.indexing import is_fancy_indexer, map_index_queries
from xarray.core.options import OPTIONS, _get_keep_attrs
from xarray.core.types import (
    Bins,
    NetcdfWriteModes,
    QuantileMethods,
    Self,
    T_ChunkDim,
    T_ChunksFreq,
    T_DataArray,
    T_DataArrayOrSet,
    ZarrWriteModes,
)
from xarray.core.utils import (
    Default,
    FilteredMapping,
    Frozen,
    FrozenMappingWarningOnValuesAccess,
    OrderedSet,
    _default,
    decode_numpy_dict_values,
    drop_dims_from_indexers,
    either_dict_or_kwargs,
    emit_user_level_warning,
    infix_dims,
    is_allowed_extension_array,
    is_dict_like,
    is_duck_array,
    is_duck_dask_array,
    is_scalar,
    maybe_wrap_array,
    parse_dims_as_set,
)
from xarray.core.variable import (
    UNSUPPORTED_EXTENSION_ARRAY_TYPES,
    IndexVariable,
    Variable,
    as_variable,
    broadcast_variables,
    calculate_dimensions,
)
from xarray.namedarray.parallelcompat import get_chunked_array_type, guess_chunkmanager
from xarray.namedarray.pycompat import array_type, is_chunked_array, to_numpy
from xarray.plot.accessor import DatasetPlotAccessor
from xarray.structure import alignment
from xarray.structure.alignment import (
    _broadcast_helper,
    _get_broadcast_dims_map_common_coords,
    align,
)
from xarray.structure.chunks import _maybe_chunk, unify_chunks
from xarray.structure.merge import (
    dataset_merge_method,
    dataset_update_method,
    merge_coordinates_without_align,
    merge_data_and_coords,
)
from xarray.util.deprecation_helpers import (
    _COMPAT_DEFAULT,
    _JOIN_DEFAULT,
    CombineKwargDefault,
    _deprecate_positional_args,
    deprecate_dims,
)

if TYPE_CHECKING:
    from dask.dataframe import DataFrame as DaskDataFrame
    from dask.delayed import Delayed
    from numpy.typing import ArrayLike

    from xarray.backends import AbstractDataStore, ZarrStore
    from xarray.backends.api import T_NetcdfEngine, T_NetcdfTypes
    from xarray.computation.rolling import DatasetCoarsen, DatasetRolling
    from xarray.computation.weighted import DatasetWeighted
    from xarray.core.dataarray import DataArray
    from xarray.core.groupby import DatasetGroupBy
    from xarray.core.resample import DatasetResample
    from xarray.core.types import (
        CFCalendar,
        CoarsenBoundaryOptions,
        CombineAttrsOptions,
        CompatOptions,
        DataVars,
        DatetimeLike,
        DatetimeUnitOptions,
        Dims,
        DsCompatible,
        ErrorOptions,
        ErrorOptionsWithWarn,
        GroupIndices,
        GroupInput,
        InterpOptions,
        JoinOptions,
        PadModeOptions,
        PadReflectOptions,
        QueryEngineOptions,
        QueryParserOptions,
        ReindexMethodOptions,
        ResampleCompatible,
        SideOptions,
        T_ChunkDimFreq,
        T_Chunks,
        T_DatasetPadConstantValues,
        T_Xarray,
    )
    from xarray.groupers import Grouper, Resampler
    from xarray.namedarray.parallelcompat import ChunkManagerEntrypoint
    from xarray.structure.merge import CoercibleMapping, CoercibleValue


# list of attributes of pd.DatetimeIndex that are ndarrays of time info
_DATETIMEINDEX_COMPONENTS = [
    "year",
    "month",
    "day",
    "hour",
    "minute",
    "second",
    "microsecond",
    "nanosecond",
    "date",
    "time",
    "dayofyear",
    "weekofyear",
    "dayofweek",
    "quarter",
]


class Dataset(
    DataWithCoords,
    DatasetAggregations,
    DatasetArithmetic,
    Mapping[Hashable, "DataArray"],
):
    """A multi-dimensional, in memory, array database.

    A dataset resembles an in-memory representation of a NetCDF file,
    and consists of variables, coordinates and attributes which
    together form a self describing dataset.

    Dataset implements the mapping interface with keys given by variable
    names and values given by DataArray objects for each variable name.

    By default, pandas indexes are created for one dimensional variables with
    name equal to their dimension (i.e., :term:`Dimension coordinate`) so those
    variables can be readily used as coordinates for label based indexing. When a
    :py:class:`~xarray.Coordinates` object is passed to ``coords``, any existing
    index(es) built from those coordinates will be added to the Dataset.

    To load data from a file or file-like object, use the `open_dataset`
    function.

    Parameters
    ----------
    data_vars : dict-like, optional
        A mapping from variable names to :py:class:`~xarray.DataArray`
        objects, :py:class:`~xarray.Variable` objects or to tuples of
        the form ``(dims, data[, attrs])`` which can be used as
        arguments to create a new ``Variable``. Each dimension must
        have the same length in all variables in which it appears.

        The following notations are accepted:

        - mapping {var name: DataArray}
        - mapping {var name: Variable}
        - mapping {var name: (dimension name, array-like)}
        - mapping {var name: (tuple of dimension names, array-like)}
        - mapping {dimension name: array-like}
          (if array-like is not a scalar it will be automatically moved to coords,
          see below)

        Each dimension must have the same length in all variables in
        which it appears.
    coords : :py:class:`~xarray.Coordinates` or dict-like, optional
        A :py:class:`~xarray.Coordinates` object or another mapping in
        similar form as the `data_vars` argument, except that each item
        is saved on the dataset as a "coordinate".
        These variables have an associated meaning: they describe
        constant/fixed/independent quantities, unlike the
        varying/measured/dependent quantities that belong in
        `variables`.

        The following notations are accepted for arbitrary mappings:

        - mapping {coord name: DataArray}
        - mapping {coord name: Variable}
        - mapping {coord name: (dimension name, array-like)}
        - mapping {coord name: (tuple of dimension names, array-like)}
        - mapping {dimension name: array-like}
          (the dimension name is implicitly set to be the same as the
          coord name)

        The last notation implies either that the coordinate value is a scalar
        or that it is a 1-dimensional array and the coord name is the same as
        the dimension name (i.e., a :term:`Dimension coordinate`). In the latter
        case, the 1-dimensional array will be assumed to give index values
        along the dimension with the same name.

        Alternatively, a :py:class:`~xarray.Coordinates` object may be used in
        order to explicitly pass indexes (e.g., a multi-index or any custom
        Xarray index) or to bypass the creation of a default index for any
        :term:`Dimension coordinate` included in that object.

    attrs : dict-like, optional
        Global attributes to save on this dataset.
        (see FAQ, :ref:`approach to metadata`)

    Examples
    --------
    In this example dataset, we will represent measurements of the temperature
    and pressure that were made under various conditions:

    * the measurements were made on four different days;
    * they were made at two separate locations, which we will represent using
      their latitude and longitude; and
    * they were made using three instrument developed by three different
      manufacturers, which we will refer to using the strings `'manufac1'`,
      `'manufac2'`, and `'manufac3'`.

    >>> np.random.seed(0)
    >>> temperature = 15 + 8 * np.random.randn(2, 3, 4)
    >>> precipitation = 10 * np.random.rand(2, 3, 4)
    >>> lon = [-99.83, -99.32]
    >>> lat = [42.25, 42.21]
    >>> instruments = ["manufac1", "manufac2", "manufac3"]
    >>> time = pd.date_range("2014-09-06", periods=4)
    >>> reference_time = pd.Timestamp("2014-09-05")

    Here, we initialize the dataset with multiple dimensions. We use the string
    `"loc"` to represent the location dimension of the data, the string
    `"instrument"` to represent the instrument manufacturer dimension, and the
    string `"time"` for the time dimension.

    >>> ds = xr.Dataset(
    ...     data_vars=dict(
    ...         temperature=(["loc", "instrument", "time"], temperature),
    ...         precipitation=(["loc", "instrument", "time"], precipitation),
    ...     ),
    ...     coords=dict(
    ...         lon=("loc", lon),
    ...         lat=("loc", lat),
    ...         instrument=instruments,
    ...         time=time,
    ...         reference_time=reference_time,
    ...     ),
    ...     attrs=dict(description="Weather related data."),
    ... )
    >>> ds
    <xarray.Dataset> Size: 552B
    Dimensions:         (loc: 2, instrument: 3, time: 4)
    Coordinates:
        lon             (loc) float64 16B -99.83 -99.32
        lat             (loc) float64 16B 42.25 42.21
      * instrument      (instrument) <U8 96B 'manufac1' 'manufac2' 'manufac3'
      * time            (time) datetime64[ns] 32B 2014-09-06 ... 2014-09-09
        reference_time  datetime64[ns] 8B 2014-09-05
    Dimensions without coordinates: loc
    Data variables:
        temperature     (loc, instrument, time) float64 192B 29.11 18.2 ... 9.063
        precipitation   (loc, instrument, time) float64 192B 4.562 5.684 ... 1.613
    Attributes:
        description:  Weather related data.

    Find out where the coldest temperature was and what values the
    other variables had:

    >>> ds.isel(ds.temperature.argmin(...))
    <xarray.Dataset> Size: 80B
    Dimensions:         ()
    Coordinates:
        lon             float64 8B -99.32
        lat             float64 8B 42.21
        instrument      <U8 32B 'manufac3'
        time            datetime64[ns] 8B 2014-09-06
        reference_time  datetime64[ns] 8B 2014-09-05
    Data variables:
        temperature     float64 8B -5.424
        precipitation   float64 8B 9.884
    Attributes:
        description:  Weather related data.

    """

    _attrs: dict[Hashable, Any] | None
    _cache: dict[str, Any]
    _coord_names: set[Hashable]
    _dims: dict[Hashable, int]
    _encoding: dict[Hashable, Any] | None
    _close: Callable[[], None] | None
    _indexes: dict[Hashable, Index]
    _variables: dict[Hashable, Variable]

    __slots__ = (
        "__weakref__",
        "_attrs",
        "_cache",
        "_close",
        "_coord_names",
        "_dims",
        "_encoding",
        "_indexes",
        "_variables",
    )

    def __init__(
        self,
        # could make a VariableArgs to use more generally, and refine these
        # categories
        data_vars: DataVars | None = None,
        coords: Mapping[Any, Any] | None = None,
        attrs: Mapping[Any, Any] | None = None,
    ) -> None:
        if data_vars is None:
            data_vars = {}
        if coords is None:
            coords = {}

        both_data_and_coords = set(data_vars) & set(coords)
        if both_data_and_coords:
            raise ValueError(
                f"variables {both_data_and_coords!r} are found in both data_vars and coords"
            )

        if isinstance(coords, Dataset):
            coords = coords._variables

        variables, coord_names, dims, indexes, _ = merge_data_and_coords(
            data_vars, coords
        )

        self._attrs = dict(attrs) if attrs else None
        self._close = None
        self._encoding = None
        self._variables = variables
        self._coord_names = coord_names
        self._dims = dims
        self._indexes = indexes

    # TODO: dirty workaround for mypy 1.5 error with inherited DatasetOpsMixin vs. Mapping
    # related to https://github.com/python/mypy/issues/9319?
    def __eq__(self, other: DsCompatible) -> Self:  # type: ignore[override]
        return super().__eq__(other)

    @classmethod
    def load_store(cls, store, decoder=None) -> Self:
        """Create a new dataset from the contents of a backends.*DataStore
        object
        """
        variables, attributes = store.load()
        if decoder:
            variables, attributes = decoder(variables, attributes)
        obj = cls(variables, attrs=attributes)
        obj.set_close(store.close)
        return obj

    @property
    def variables(self) -> Frozen[Hashable, Variable]:
        """Low level interface to Dataset contents as dict of Variable objects.

        This ordered dictionary is frozen to prevent mutation that could
        violate Dataset invariants. It contains all variable objects
        constituting the Dataset, including both data variables and
        coordinates.
        """
        return Frozen(self._variables)

    @property
    def attrs(self) -> dict[Any, Any]:
        """Dictionary of global attributes on this dataset"""
        if self._attrs is None:
            self._attrs = {}
        return self._attrs

    @attrs.setter
    def attrs(self, value: Mapping[Any, Any]) -> None:
        self._attrs = dict(value) if value else None

    @property
    def encoding(self) -> dict[Any, Any]:
        """Dictionary of global encoding attributes on this dataset"""
        if self._encoding is None:
            self._encoding = {}
        return self._encoding

    @encoding.setter
    def encoding(self, value: Mapping[Any, Any]) -> None:
        self._encoding = dict(value)

    def reset_encoding(self) -> Self:
        warnings.warn(
            "reset_encoding is deprecated since 2023.11, use `drop_encoding` instead",
            stacklevel=2,
        )
        return self.drop_encoding()

    def drop_encoding(self) -> Self:
        """Return a new Dataset without encoding on the dataset or any of its
        variables/coords."""
        variables = {k: v.drop_encoding() for k, v in self.variables.items()}
        return self._replace(variables=variables, encoding={})

    @property
    def dims(self) -> Frozen[Hashable, int]:
        """Mapping from dimension names to lengths.

        Cannot be modified directly, but is updated when adding new variables.

        Note that type of this object differs from `DataArray.dims`.
        See `Dataset.sizes` and `DataArray.sizes` for consistently named
        properties. This property will be changed to return a type more consistent with
        `DataArray.dims` in the future, i.e. a set of dimension names.

        See Also
        --------
        Dataset.sizes
        DataArray.dims
        """
        return FrozenMappingWarningOnValuesAccess(self._dims)

    @property
    def sizes(self) -> Frozen[Hashable, int]:
        """Mapping from dimension names to lengths.

        Cannot be modified directly, but is updated when adding new variables.

        This is an alias for `Dataset.dims` provided for the benefit of
        consistency with `DataArray.sizes`.

        See Also
        --------
        DataArray.sizes
        """
        return Frozen(self._dims)

    @property
    def dtypes(self) -> Frozen[Hashable, np.dtype]:
        """Mapping from data variable names to dtypes.

        Cannot be modified directly, but is updated when adding new variables.

        See Also
        --------
        DataArray.dtype
        """
        return Frozen(
            {
                n: v.dtype
                for n, v in self._variables.items()
                if n not in self._coord_names
            }
        )

    def load(self, **kwargs) -> Self:
        """Trigger loading data into memory and return this dataset.

        Data will be computed and/or loaded from disk or a remote source.

        Unlike ``.compute``, the original dataset is modified and returned.

        Normally, it should not be necessary to call this method in user code,
        because all xarray functions should either work on deferred data or
        load data automatically. However, this method can be necessary when
        working with many file objects on disk.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.compute``.

        Returns
        -------
        object : Dataset
            Same object but with lazy data variables and coordinates as in-memory arrays.

        See Also
        --------
        dask.compute
        Dataset.compute
        Dataset.load_async
        DataArray.load
        Variable.load
        """
        # access .data to coerce everything to numpy or dask arrays
        chunked_data = {
            k: v._data for k, v in self.variables.items() if is_chunked_array(v._data)
        }
        if chunked_data:
            chunkmanager = get_chunked_array_type(*chunked_data.values())

            # evaluate all the chunked arrays simultaneously
            evaluated_data: tuple[np.ndarray[Any, Any], ...] = chunkmanager.compute(
                *chunked_data.values(), **kwargs
            )

            for k, data in zip(chunked_data, evaluated_data, strict=False):
                self.variables[k].data = data

        # load everything else sequentially
        [v.load() for k, v in self.variables.items() if k not in chunked_data]

        return self

    async def load_async(self, **kwargs) -> Self:
        """Trigger and await asynchronous loading of data into memory and return this dataset.

        Data will be computed and/or loaded from disk or a remote source.

        Unlike ``.compute``, the original dataset is modified and returned.

        Only works when opening data lazily from IO storage backends which support lazy asynchronous loading.
        Otherwise will raise a NotImplementedError.

        Note users are expected to limit concurrency themselves - xarray does not internally limit concurrency in any way.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.compute``.

        Returns
        -------
        object : Dataset
            Same object but with lazy data variables and coordinates as in-memory arrays.

        See Also
        --------
        dask.compute
        Dataset.compute
        Dataset.load
        DataArray.load_async
        Variable.load_async
        """
        # TODO refactor this to pull out the common chunked_data codepath

        # this blocks on chunked arrays but not on lazily indexed arrays

        # access .data to coerce everything to numpy or dask arrays
        chunked_data = {
            k: v._data for k, v in self.variables.items() if is_chunked_array(v._data)
        }
        if chunked_data:
            chunkmanager = get_chunked_array_type(*chunked_data.values())

            # evaluate all the chunked arrays simultaneously
            evaluated_data: tuple[np.ndarray[Any, Any], ...] = chunkmanager.compute(
                *chunked_data.values(), **kwargs
            )

            for k, data in zip(chunked_data, evaluated_data, strict=False):
                self.variables[k].data = data

        # load everything else concurrently
        coros = [
            v.load_async() for k, v in self.variables.items() if k not in chunked_data
        ]
        await asyncio.gather(*coros)

        return self

    def __dask_tokenize__(self) -> object:
        from dask.base import normalize_token

        return normalize_token(
            (type(self), self._variables, self._coord_names, self._attrs or None)
        )

    def __dask_graph__(self):
        graphs = {k: v.__dask_graph__() for k, v in self.variables.items()}
        graphs = {k: v for k, v in graphs.items() if v is not None}
        if not graphs:
            return None
        else:
            try:
                from dask.highlevelgraph import HighLevelGraph

                return HighLevelGraph.merge(*graphs.values())
            except ImportError:
                from dask import sharedict

                return sharedict.merge(*graphs.values())

    def __dask_keys__(self):
        import dask

        return [
            v.__dask_keys__()
            for v in self.variables.values()
            if dask.is_dask_collection(v)
        ]

    def __dask_layers__(self):
        import dask

        return sum(
            (
                v.__dask_layers__()
                for v in self.variables.values()
                if dask.is_dask_collection(v)
            ),
            (),
        )

    @property
    def __dask_optimize__(self):
        import dask.array as da

        return da.Array.__dask_optimize__

    @property
    def __dask_scheduler__(self):
        import dask.array as da

        return da.Array.__dask_scheduler__

    def __dask_postcompute__(self):
        return self._dask_postcompute, ()

    def __dask_postpersist__(self):
        return self._dask_postpersist, ()

    def _dask_postcompute(self, results: Iterable[Variable]) -> Self:
        import dask

        variables = {}
        results_iter = iter(results)

        for k, v in self._variables.items():
            if dask.is_dask_collection(v):
                rebuild, args = v.__dask_postcompute__()
                v = rebuild(next(results_iter), *args)
            variables[k] = v

        return type(self)._construct_direct(
            variables,
            self._coord_names,
            self._dims,
            self._attrs,
            self._indexes,
            self._encoding,
            self._close,
        )

    def _dask_postpersist(
        self, dsk: Mapping, *, rename: Mapping[str, str] | None = None
    ) -> Self:
        from dask import is_dask_collection
        from dask.highlevelgraph import HighLevelGraph
        from dask.optimization import cull

        variables = {}

        for k, v in self._variables.items():
            if not is_dask_collection(v):
                variables[k] = v
                continue

            if isinstance(dsk, HighLevelGraph):
                # dask >= 2021.3
                # __dask_postpersist__() was called by dask.highlevelgraph.
                # Don't use dsk.cull(), as we need to prevent partial layers:
                # https://github.com/dask/dask/issues/7137
                layers = v.__dask_layers__()
                if rename:
                    layers = [rename.get(k, k) for k in layers]
                dsk2 = dsk.cull_layers(layers)
            elif rename:  # pragma: nocover
                # At the moment of writing, this is only for forward compatibility.
                # replace_name_in_key requires dask >= 2021.3.
                from dask.base import flatten, replace_name_in_key

                keys = [
                    replace_name_in_key(k, rename) for k in flatten(v.__dask_keys__())
                ]
                dsk2, _ = cull(dsk, keys)
            else:
                # __dask_postpersist__() was called by dask.optimize or dask.persist
                dsk2, _ = cull(dsk, v.__dask_keys__())

            rebuild, args = v.__dask_postpersist__()
            # rename was added in dask 2021.3
            kwargs = {"rename": rename} if rename else {}
            variables[k] = rebuild(dsk2, *args, **kwargs)

        return type(self)._construct_direct(
            variables,
            self._coord_names,
            self._dims,
            self._attrs,
            self._indexes,
            self._encoding,
            self._close,
        )

    def compute(self, **kwargs) -> Self:
        """Trigger loading data into memory and return a new dataset.

        Data will be computed and/or loaded from disk or a remote source.

        Unlike ``.load``, the original dataset is left unaltered.

        Normally, it should not be necessary to call this method in user code,
        because all xarray functions should either work on deferred data or
        load data automatically. However, this method can be necessary when
        working with many file objects on disk.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.compute``.

        Returns
        -------
        object : Dataset
            New object with lazy data variables and coordinates as in-memory arrays.

        See Also
        --------
        dask.compute
        Dataset.load
        Dataset.load_async
        DataArray.compute
        Variable.compute
        """
        new = self.copy(deep=False)
        return new.load(**kwargs)

    def _persist_inplace(self, **kwargs) -> Self:
        """Persist all chunked arrays in memory."""
        # access .data to coerce everything to numpy or dask arrays
        lazy_data = {
            k: v._data for k, v in self.variables.items() if is_chunked_array(v._data)
        }
        if lazy_data:
            chunkmanager = get_chunked_array_type(*lazy_data.values())

            # evaluate all the dask arrays simultaneously
            evaluated_data = chunkmanager.persist(*lazy_data.values(), **kwargs)

            for k, data in zip(lazy_data, evaluated_data, strict=False):
                self.variables[k].data = data

        return self

    def persist(self, **kwargs) -> Self:
        """Trigger computation, keeping data as chunked arrays.

        This operation can be used to trigger computation on underlying dask
        arrays, similar to ``.compute()`` or ``.load()``.  However this
        operation keeps the data as dask arrays. This is particularly useful
        when using the dask.distributed scheduler and you want to load a large
        amount of data into distributed memory.
        Like compute (but unlike load), the original dataset is left unaltered.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.persist``.

        Returns
        -------
        object : Dataset
            New object with all dask-backed coordinates and data variables as persisted dask arrays.

        See Also
        --------
        dask.persist
        """
        new = self.copy(deep=False)
        return new._persist_inplace(**kwargs)

    @classmethod
    def _construct_direct(
        cls,
        variables: dict[Any, Variable],
        coord_names: set[Hashable],
        dims: dict[Any, int] | None = None,
        attrs: dict | None = None,
        indexes: dict[Any, Index] | None = None,
        encoding: dict | None = None,
        close: Callable[[], None] | None = None,
    ) -> Self:
        """Shortcut around __init__ for internal use when we want to skip
        costly validation
        """
        if dims is None:
            dims = calculate_dimensions(variables)
        if indexes is None:
            indexes = {}
        obj = object.__new__(cls)
        obj._variables = variables
        obj._coord_names = coord_names
        obj._dims = dims
        obj._indexes = indexes
        obj._attrs = attrs
        obj._close = close
        obj._encoding = encoding
        return obj

    def _replace(
        self,
        variables: dict[Hashable, Variable] | None = None,
        coord_names: set[Hashable] | None = None,
        dims: dict[Any, int] | None = None,
        attrs: dict[Hashable, Any] | Default | None = _default,
        indexes: dict[Hashable, Index] | None = None,
        encoding: dict | Default | None = _default,
        inplace: bool = False,
    ) -> Self:
        """Fastpath constructor for internal use.

        Returns an object with optionally with replaced attributes.

        Explicitly passed arguments are *not* copied when placed on the new
        dataset. It is up to the caller to ensure that they have the right type
        and are not used elsewhere.
        """
        if inplace:
            if variables is not None:
                self._variables = variables
            if coord_names is not None:
                self._coord_names = coord_names
            if dims is not None:
                self._dims = dims
            if attrs is not _default:
                self._attrs = attrs
            if indexes is not None:
                self._indexes = indexes
            if encoding is not _default:
                self._encoding = encoding
            obj = self
        else:
            if variables is None:
                variables = self._variables.copy()
            if coord_names is None:
                coord_names = self._coord_names.copy()
            if dims is None:
                dims = self._dims.copy()
            if attrs is _default:
                attrs = copy.copy(self._attrs)
            if indexes is None:
                indexes = self._indexes.copy()
            if encoding is _default:
                encoding = copy.copy(self._encoding)
            obj = self._construct_direct(
                variables, coord_names, dims, attrs, indexes, encoding
            )
        return obj

    def _replace_with_new_dims(
        self,
        variables: dict[Hashable, Variable],
        coord_names: set | None = None,
        attrs: dict[Hashable, Any] | Default | None = _default,
        indexes: dict[Hashable, Index] | None = None,
        inplace: bool = False,
    ) -> Self:
        """Replace variables with recalculated dimensions."""
        dims = calculate_dimensions(variables)
        return self._replace(
            variables, coord_names, dims, attrs, indexes, inplace=inplace
        )

    def _replace_vars_and_dims(
        self,
        variables: dict[Hashable, Variable],
        coord_names: set | None = None,
        dims: dict[Hashable, int] | None = None,
        attrs: dict[Hashable, Any] | Default | None = _default,
        inplace: bool = False,
    ) -> Self:
        """Deprecated version of _replace_with_new_dims().

        Unlike _replace_with_new_dims(), this method always recalculates
        indexes from variables.
        """
        if dims is None:
            dims = calculate_dimensions(variables)
        return self._replace(
            variables, coord_names, dims, attrs, indexes=None, inplace=inplace
        )

    def _overwrite_indexes(
        self,
        indexes: Mapping[Hashable, Index],
        variables: Mapping[Hashable, Variable] | None = None,
        drop_variables: list[Hashable] | None = None,
        drop_indexes: list[Hashable] | None = None,
        rename_dims: Mapping[Hashable, Hashable] | None = None,
    ) -> Self:
        """Maybe replace indexes.

        This function may do a lot more depending on index query
        results.

        """
        if not indexes:
            return self

        if variables is None:
            variables = {}
        if drop_variables is None:
            drop_variables = []
        if drop_indexes is None:
            drop_indexes = []

        new_variables = self._variables.copy()
        new_coord_names = self._coord_names.copy()
        new_indexes = dict(self._indexes)

        index_variables = {}
        no_index_variables = {}
        for name, var in variables.items():
            old_var = self._variables.get(name)
            if old_var is not None:
                var.attrs.update(old_var.attrs)
                var.encoding.update(old_var.encoding)
            if name in indexes:
                index_variables[name] = var
            else:
                no_index_variables[name] = var

        for name in indexes:
            new_indexes[name] = indexes[name]

        for name, var in index_variables.items():
            new_coord_names.add(name)
            new_variables[name] = var

        # append no-index variables at the end
        for k in no_index_variables:
            new_variables.pop(k)
        new_variables.update(no_index_variables)

        for name in drop_indexes:
            new_indexes.pop(name)

        for name in drop_variables:
            new_variables.pop(name)
            new_indexes.pop(name, None)
            new_coord_names.remove(name)

        replaced = self._replace(
            variables=new_variables, coord_names=new_coord_names, indexes=new_indexes
        )

        if rename_dims:
            # skip rename indexes: they should already have the right name(s)
            dims = replaced._rename_dims(rename_dims)
            new_variables, new_coord_names = replaced._rename_vars({}, rename_dims)
            return replaced._replace(
                variables=new_variables, coord_names=new_coord_names, dims=dims
            )
        else:
            return replaced

    def copy(self, deep: bool = False, data: DataVars | None = None) -> Self:
        """Returns a copy of this dataset.

        If `deep=True`, a deep copy is made of each of the component variables.
        Otherwise, a shallow copy of each of the component variable is made, so
        that the underlying memory region of the new dataset is the same as in
        the original dataset.

        Use `data` to create a new object with the same structure as
        original but entirely new data.

        Parameters
        ----------
        deep : bool, default: False
            Whether each component variable is loaded into memory and copied onto
            the new object. Default is False.
        data : dict-like or None, optional
            Data to use in the new object. Each item in `data` must have same
            shape as corresponding data variable in original. When `data` is
            used, `deep` is ignored for the data variables and only used for
            coords.

        Returns
        -------
        object : Dataset
            New object with dimensions, attributes, coordinates, name, encoding,
            and optionally data copied from original.

        Examples
        --------
        Shallow copy versus deep copy

        >>> da = xr.DataArray(np.random.randn(2, 3))
        >>> ds = xr.Dataset(
        ...     {"foo": da, "bar": ("x", [-1, 2])},
        ...     coords={"x": ["one", "two"]},
        ... )
        >>> ds.copy()
        <xarray.Dataset> Size: 88B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Coordinates:
          * x        (x) <U3 24B 'one' 'two'
        Dimensions without coordinates: dim_0, dim_1
        Data variables:
            foo      (dim_0, dim_1) float64 48B 1.764 0.4002 0.9787 2.241 1.868 -0.9773
            bar      (x) int64 16B -1 2

        >>> ds_0 = ds.copy(deep=False)
        >>> ds_0["foo"][0, 0] = 7
        >>> ds_0
        <xarray.Dataset> Size: 88B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Coordinates:
          * x        (x) <U3 24B 'one' 'two'
        Dimensions without coordinates: dim_0, dim_1
        Data variables:
            foo      (dim_0, dim_1) float64 48B 7.0 0.4002 0.9787 2.241 1.868 -0.9773
            bar      (x) int64 16B -1 2

        >>> ds
        <xarray.Dataset> Size: 88B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Coordinates:
          * x        (x) <U3 24B 'one' 'two'
        Dimensions without coordinates: dim_0, dim_1
        Data variables:
            foo      (dim_0, dim_1) float64 48B 7.0 0.4002 0.9787 2.241 1.868 -0.9773
            bar      (x) int64 16B -1 2

        Changing the data using the ``data`` argument maintains the
        structure of the original object, but with the new data. Original
        object is unaffected.

        >>> ds.copy(data={"foo": np.arange(6).reshape(2, 3), "bar": ["a", "b"]})
        <xarray.Dataset> Size: 80B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Coordinates:
          * x        (x) <U3 24B 'one' 'two'
        Dimensions without coordinates: dim_0, dim_1
        Data variables:
            foo      (dim_0, dim_1) int64 48B 0 1 2 3 4 5
            bar      (x) <U1 8B 'a' 'b'

        >>> ds
        <xarray.Dataset> Size: 88B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Coordinates:
          * x        (x) <U3 24B 'one' 'two'
        Dimensions without coordinates: dim_0, dim_1
        Data variables:
            foo      (dim_0, dim_1) float64 48B 7.0 0.4002 0.9787 2.241 1.868 -0.9773
            bar      (x) int64 16B -1 2

        See Also
        --------
        pandas.DataFrame.copy
        """
        return self._copy(deep=deep, data=data)

    def _copy(
        self,
        deep: bool = False,
        data: DataVars | None = None,
        memo: dict[int, Any] | None = None,
    ) -> Self:
        if data is None:
            data = {}
        elif not utils.is_dict_like(data):
            raise ValueError("Data must be dict-like")

        if data:
            var_keys = set(self.data_vars.keys())
            data_keys = set(data.keys())
            keys_not_in_vars = data_keys - var_keys
            if keys_not_in_vars:
                raise ValueError(
                    "Data must only contain variables in original "
                    f"dataset. Extra variables: {keys_not_in_vars}"
                )
            keys_missing_from_data = var_keys - data_keys
            if keys_missing_from_data:
                raise ValueError(
                    "Data must contain all variables in original "
                    f"dataset. Data is missing {keys_missing_from_data}"
                )

        indexes, index_vars = self.xindexes.copy_indexes(deep=deep)

        variables = {}
        for k, v in self._variables.items():
            if k in index_vars:
                variables[k] = index_vars[k]
            else:
                variables[k] = v._copy(deep=deep, data=data.get(k), memo=memo)

        attrs = copy.deepcopy(self._attrs, memo) if deep else copy.copy(self._attrs)
        encoding = (
            copy.deepcopy(self._encoding, memo) if deep else copy.copy(self._encoding)
        )

        return self._replace(variables, indexes=indexes, attrs=attrs, encoding=encoding)

    def __copy__(self) -> Self:
        return self._copy(deep=False)

    def __deepcopy__(self, memo: dict[int, Any] | None = None) -> Self:
        return self._copy(deep=True, memo=memo)

    def as_numpy(self) -> Self:
        """
        Coerces wrapped data and coordinates into numpy arrays, returning a Dataset.

        See also
        --------
        DataArray.as_numpy
        DataArray.to_numpy : Returns only the data as a numpy.ndarray object.
        """
        numpy_variables = {k: v.as_numpy() for k, v in self.variables.items()}
        return self._replace(variables=numpy_variables)

    def _copy_listed(self, names: Iterable[Hashable]) -> Self:
        """Create a new Dataset with the listed variables from this dataset and
        the all relevant coordinates. Skips all validation.
        """
        variables: dict[Hashable, Variable] = {}
        coord_names = set()
        indexes: dict[Hashable, Index] = {}

        for name in names:
            try:
                variables[name] = self._variables[name]
            except KeyError:
                ref_name, var_name, var = _get_virtual_variable(
                    self._variables, name, self.sizes
                )
                variables[var_name] = var
                if ref_name in self._coord_names or ref_name in self.dims:
                    coord_names.add(var_name)
                if (var_name,) == var.dims:
                    index, index_vars = create_default_index_implicit(var, names)
                    indexes.update(dict.fromkeys(index_vars, index))
                    variables.update(index_vars)
                    coord_names.update(index_vars)

        needed_dims: OrderedSet[Hashable] = OrderedSet()
        for v in variables.values():
            needed_dims.update(v.dims)

        dims = {k: self.sizes[k] for k in needed_dims}

        # preserves ordering of coordinates
        for k in self._variables:
            if k not in self._coord_names:
                continue

            if set(self.variables[k].dims) <= needed_dims:
                variables[k] = self._variables[k]
                coord_names.add(k)

        indexes.update(filter_indexes_from_coords(self._indexes, coord_names))

        return self._replace(variables, coord_names, dims, indexes=indexes)

    def _construct_dataarray(self, name: Hashable) -> DataArray:
        """Construct a DataArray by indexing this dataset"""
        from xarray.core.dataarray import DataArray

        try:
            variable = self._variables[name]
        except KeyError:
            _, name, variable = _get_virtual_variable(self._variables, name, self.sizes)

        needed_dims = set(variable.dims)

        coords: dict[Hashable, Variable] = {}
        # preserve ordering
        for k in self._variables:
            if k in self._indexes:
                add_coord = self._indexes[k].should_add_coord_to_array(
                    k, self._variables[k], needed_dims
                )
            else:
                var_dims = set(self._variables[k].dims)
                add_coord = k in self._coord_names and var_dims <= needed_dims

            if add_coord:
                coords[k] = self._variables[k]

        indexes = filter_indexes_from_coords(self._indexes, set(coords))

        return DataArray(variable, coords, name=name, indexes=indexes, fastpath=True)

    @property
    def _attr_sources(self) -> Iterable[Mapping[Hashable, Any]]:
        """Places to look-up items for attribute-style access"""
        yield from self._item_sources
        yield self.attrs

    @property
    def _item_sources(self) -> Iterable[Mapping[Hashable, Any]]:
        """Places to look-up items for key-completion"""
        yield self.data_vars
        yield FilteredMapping(keys=self._coord_names, mapping=self.coords)

        # virtual coordinates
        yield FilteredMapping(keys=self.sizes, mapping=self)

    def __contains__(self, key: object) -> bool:
        """The 'in' operator will return true or false depending on whether
        'key' is an array in the dataset or not.
        """
        return key in self._variables

    def __len__(self) -> int:
        return len(self.data_vars)

    def __bool__(self) -> bool:
        return bool(self.data_vars)

    def __iter__(self) -> Iterator[Hashable]:
        return iter(self.data_vars)

    if TYPE_CHECKING:
        # needed because __getattr__ is returning Any and otherwise
        # this class counts as part of the SupportsArray Protocol
        __array__ = None  # type: ignore[var-annotated,unused-ignore]

    else:

        def __array__(self, dtype=None, copy=None):
            raise TypeError(
                "cannot directly convert an xarray.Dataset into a "
                "numpy array. Instead, create an xarray.DataArray "
                "first, either with indexing on the Dataset or by "
                "invoking the `to_dataarray()` method."
            )

    @property
    def nbytes(self) -> int:
        """
        Total bytes consumed by the data arrays of all variables in this dataset.

        If the backend array for any variable does not include ``nbytes``, estimates
        the total bytes for that array based on the ``size`` and ``dtype``.
        """
        return sum(v.nbytes for v in self.variables.values())

    @property
    def loc(self) -> _LocIndexer[Self]:
        """Attribute for location based indexing. Only supports __getitem__,
        and only when the key is a dict of the form {dim: labels}.
        """
        return _LocIndexer(self)

    @overload
    def __getitem__(self, key: Hashable) -> DataArray: ...

    # Mapping is Iterable
    @overload
    def __getitem__(self, key: Iterable[Hashable]) -> Self: ...

    def __getitem__(
        self, key: Mapping[Any, Any] | Hashable | Iterable[Hashable]
    ) -> Self | DataArray:
        """Access variables or coordinates of this dataset as a
        :py:class:`~xarray.DataArray` or a subset of variables or a indexed dataset.

        Indexing with a list of names will return a new ``Dataset`` object.
        """
        from xarray.core.formatting import shorten_list_repr

        if utils.is_dict_like(key):
            return self.isel(**key)
        if utils.hashable(key):
            try:
                return self._construct_dataarray(key)
            except KeyError as e:
                message = f"No variable named {key!r}."

                best_guess = utils.did_you_mean(key, self.variables.keys())
                if best_guess:
                    message += f" {best_guess}"
                else:
                    message += f" Variables on the dataset include {shorten_list_repr(list(self.variables.keys()), max_items=10)}"

                # If someone attempts `ds['foo' , 'bar']` instead of `ds[['foo', 'bar']]`
                if isinstance(key, tuple):
                    message += f"\nHint: use a list to select multiple variables, for example `ds[{list(key)}]`"
                raise KeyError(message) from e

        if utils.iterable_of_hashable(key):
            return self._copy_listed(key)
        raise ValueError(f"Unsupported key-type {type(key)}")

    def __setitem__(
        self, key: Hashable | Iterable[Hashable] | Mapping, value: Any
    ) -> None:
        """Add an array to this dataset.
        Multiple arrays can be added at the same time, in which case each of
        the following operations is applied to the respective value.

        If key is dict-like, update all variables in the dataset
        one by one with the given value at the given location.
        If the given value is also a dataset, select corresponding variables
        in the given value and in the dataset to be changed.

        If value is a `
        from .dataarray import DataArray`, call its `select_vars()` method, rename it
        to `key` and merge the contents of the resulting dataset into this
        dataset.

        If value is a `Variable` object (or tuple of form
        ``(dims, data[, attrs])``), add it to this dataset as a new
        variable.
        """
        from xarray.core.dataarray import DataArray

        if utils.is_dict_like(key):
            # check for consistency and convert value to dataset
            value = self._setitem_check(key, value)
            # loop over dataset variables and set new values
            processed = []
            for name, var in self.items():
                try:
                    var[key] = value[name]
                    processed.append(name)
                except Exception as e:
                    if processed:
                        raise RuntimeError(
                            "An error occurred while setting values of the"
                            f" variable '{name}'. The following variables have"
                            f" been successfully updated:\n{processed}"
                        ) from e
                    else:
                        raise e

        elif utils.hashable(key):
            if isinstance(value, Dataset):
                raise TypeError(
                    "Cannot assign a Dataset to a single key - only a DataArray or Variable "
                    "object can be stored under a single key."
                )
            self.update({key: value})

        elif utils.iterable_of_hashable(key):
            keylist = list(key)
            if len(keylist) == 0:
                raise ValueError("Empty list of variables to be set")
            if len(keylist) == 1:
                self.update({keylist[0]: value})
            else:
                if len(keylist) != len(value):
                    raise ValueError(
                        f"Different lengths of variables to be set "
                        f"({len(keylist)}) and data used as input for "
                        f"setting ({len(value)})"
                    )
                if isinstance(value, Dataset):
                    self.update(
                        dict(zip(keylist, value.data_vars.values(), strict=True))
                    )
                elif isinstance(value, DataArray):
                    raise ValueError("Cannot assign single DataArray to multiple keys")
                else:
                    self.update(dict(zip(keylist, value, strict=True)))

        else:
            raise ValueError(f"Unsupported key-type {type(key)}")

    def _setitem_check(self, key, value):
        """Consistency check for __setitem__

        When assigning values to a subset of a Dataset, do consistency check beforehand
        to avoid leaving the dataset in a partially updated state when an error occurs.
        """
        from xarray.core.dataarray import DataArray

        if isinstance(value, Dataset):
            missing_vars = [
                name for name in value.data_vars if name not in self.data_vars
            ]
            if missing_vars:
                raise ValueError(
                    f"Variables {missing_vars} in new values"
                    f" not available in original dataset:\n{self}"
                )
        elif not any(isinstance(value, t) for t in [DataArray, Number, str]):
            raise TypeError(
                "Dataset assignment only accepts DataArrays, Datasets, and scalars."
            )

        new_value = Dataset()
        for name, var in self.items():
            # test indexing
            try:
                var_k = var[key]
            except Exception as e:
                raise ValueError(
                    f"Variable '{name}': indexer {key} not available"
                ) from e

            if isinstance(value, Dataset):
                val = value[name]
            else:
                val = value

            if isinstance(val, DataArray):
                # check consistency of dimensions
                for dim in val.dims:
                    if dim not in var_k.dims:
                        raise KeyError(
                            f"Variable '{name}': dimension '{dim}' appears in new values "
                            f"but not in the indexed original data"
                        )
                dims = tuple(dim for dim in var_k.dims if dim in val.dims)
                if dims != val.dims:
                    raise ValueError(
                        f"Variable '{name}': dimension order differs between"
                        f" original and new data:\n{dims}\nvs.\n{val.dims}"
                    )
            else:
                val = np.array(val)

            # type conversion
            new_value[name] = duck_array_ops.astype(val, dtype=var_k.dtype, copy=False)

        # check consistency of dimension sizes and dimension coordinates
        if isinstance(value, DataArray | Dataset):
            align(self[key], value, join="exact", copy=False)

        return new_value

    def __delitem__(self, key: Hashable) -> None:
        """Remove a variable from this dataset."""
        assert_no_index_corrupted(self.xindexes, {key})

        if key in self._indexes:
            del self._indexes[key]
        del self._variables[key]
        self._coord_names.discard(key)
        self._dims = calculate_dimensions(self._variables)

    # mutable objects should not be hashable
    # https://github.com/python/mypy/issues/4266
    __hash__ = None  # type: ignore[assignment]

    def _all_compat(self, other: Self, compat_str: str) -> bool:
        """Helper function for equals and identical"""

        # some stores (e.g., scipy) do not seem to preserve order, so don't
        # require matching order for equality
        def compat(x: Variable, y: Variable) -> bool:
            return getattr(x, compat_str)(y)

        return self._coord_names == other._coord_names and utils.dict_equiv(
            self._variables, other._variables, compat=compat
        )

    def broadcast_equals(self, other: Self) -> bool:
        """Two Datasets are broadcast equal if they are equal after
        broadcasting all variables against each other.

        For example, variables that are scalar in one dataset but non-scalar in
        the other dataset can still be broadcast equal if the the non-scalar
        variable is a constant.

        Examples
        --------

        # 2D array with shape (1, 3)

        >>> data = np.array([[1, 2, 3]])
        >>> a = xr.Dataset(
        ...     {"variable_name": (("space", "time"), data)},
        ...     coords={"space": [0], "time": [0, 1, 2]},
        ... )
        >>> a
        <xarray.Dataset> Size: 56B
        Dimensions:        (space: 1, time: 3)
        Coordinates:
          * space          (space) int64 8B 0
          * time           (time) int64 24B 0 1 2
        Data variables:
            variable_name  (space, time) int64 24B 1 2 3

        # 2D array with shape (3, 1)

        >>> data = np.array([[1], [2], [3]])
        >>> b = xr.Dataset(
        ...     {"variable_name": (("time", "space"), data)},
        ...     coords={"time": [0, 1, 2], "space": [0]},
        ... )
        >>> b
        <xarray.Dataset> Size: 56B
        Dimensions:        (time: 3, space: 1)
        Coordinates:
          * time           (time) int64 24B 0 1 2
          * space          (space) int64 8B 0
        Data variables:
            variable_name  (time, space) int64 24B 1 2 3

        .equals returns True if two Datasets have the same values, dimensions, and coordinates. .broadcast_equals returns True if the
        results of broadcasting two Datasets against each other have the same values, dimensions, and coordinates.

        >>> a.equals(b)
        False

        >>> a.broadcast_equals(b)
        True

        >>> a2, b2 = xr.broadcast(a, b)
        >>> a2.equals(b2)
        True

        See Also
        --------
        Dataset.equals
        Dataset.identical
        Dataset.broadcast
        """
        try:
            return self._all_compat(other, "broadcast_equals")
        except (TypeError, AttributeError):
            return False

    def equals(self, other: Self) -> bool:
        """Two Datasets are equal if they have matching variables and
        coordinates, all of which are equal.

        Datasets can still be equal (like pandas objects) if they have NaN
        values in the same locations.

        This method is necessary because `v1 == v2` for ``Dataset``
        does element-wise comparisons (like numpy.ndarrays).

        Examples
        --------

        # 2D array with shape (1, 3)

        >>> data = np.array([[1, 2, 3]])
        >>> dataset1 = xr.Dataset(
        ...     {"variable_name": (("space", "time"), data)},
        ...     coords={"space": [0], "time": [0, 1, 2]},
        ... )
        >>> dataset1
        <xarray.Dataset> Size: 56B
        Dimensions:        (space: 1, time: 3)
        Coordinates:
          * space          (space) int64 8B 0
          * time           (time) int64 24B 0 1 2
        Data variables:
            variable_name  (space, time) int64 24B 1 2 3

        # 2D array with shape (3, 1)

        >>> data = np.array([[1], [2], [3]])
        >>> dataset2 = xr.Dataset(
        ...     {"variable_name": (("time", "space"), data)},
        ...     coords={"time": [0, 1, 2], "space": [0]},
        ... )
        >>> dataset2
        <xarray.Dataset> Size: 56B
        Dimensions:        (time: 3, space: 1)
        Coordinates:
          * time           (time) int64 24B 0 1 2
          * space          (space) int64 8B 0
        Data variables:
            variable_name  (time, space) int64 24B 1 2 3
        >>> dataset1.equals(dataset2)
        False

        >>> dataset1.broadcast_equals(dataset2)
        True

        .equals returns True if two Datasets have the same values, dimensions, and coordinates. .broadcast_equals returns True if the
        results of broadcasting two Datasets against each other have the same values, dimensions, and coordinates.

        Similar for missing values too:

        >>> ds1 = xr.Dataset(
        ...     {
        ...         "temperature": (["x", "y"], [[1, np.nan], [3, 4]]),
        ...     },
        ...     coords={"x": [0, 1], "y": [0, 1]},
        ... )

        >>> ds2 = xr.Dataset(
        ...     {
        ...         "temperature": (["x", "y"], [[1, np.nan], [3, 4]]),
        ...     },
        ...     coords={"x": [0, 1], "y": [0, 1]},
        ... )
        >>> ds1.equals(ds2)
        True

        See Also
        --------
        Dataset.broadcast_equals
        Dataset.identical
        """
        try:
            return self._all_compat(other, "equals")
        except (TypeError, AttributeError):
            return False

    def identical(self, other: Self) -> bool:
        """Like equals, but also checks all dataset attributes and the
        attributes on all variables and coordinates.

        Example
        -------

        >>> a = xr.Dataset(
        ...     {"Width": ("X", [1, 2, 3])},
        ...     coords={"X": [1, 2, 3]},
        ...     attrs={"units": "m"},
        ... )
        >>> b = xr.Dataset(
        ...     {"Width": ("X", [1, 2, 3])},
        ...     coords={"X": [1, 2, 3]},
        ...     attrs={"units": "m"},
        ... )
        >>> c = xr.Dataset(
        ...     {"Width": ("X", [1, 2, 3])},
        ...     coords={"X": [1, 2, 3]},
        ...     attrs={"units": "ft"},
        ... )
        >>> a
        <xarray.Dataset> Size: 48B
        Dimensions:  (X: 3)
        Coordinates:
          * X        (X) int64 24B 1 2 3
        Data variables:
            Width    (X) int64 24B 1 2 3
        Attributes:
            units:    m

        >>> b
        <xarray.Dataset> Size: 48B
        Dimensions:  (X: 3)
        Coordinates:
          * X        (X) int64 24B 1 2 3
        Data variables:
            Width    (X) int64 24B 1 2 3
        Attributes:
            units:    m

        >>> c
        <xarray.Dataset> Size: 48B
        Dimensions:  (X: 3)
        Coordinates:
          * X        (X) int64 24B 1 2 3
        Data variables:
            Width    (X) int64 24B 1 2 3
        Attributes:
            units:    ft

        >>> a.equals(b)
        True

        >>> a.identical(b)
        True

        >>> a.equals(c)
        True

        >>> a.identical(c)
        False

        See Also
        --------
        Dataset.broadcast_equals
        Dataset.equals
        """
        try:
            return utils.dict_equiv(self.attrs, other.attrs) and self._all_compat(
                other, "identical"
            )
        except (TypeError, AttributeError):
            return False

    @property
    def indexes(self) -> Indexes[pd.Index]:
        """Mapping of pandas.Index objects used for label based indexing.

        Raises an error if this Dataset has indexes that cannot be coerced
        to pandas.Index objects.

        See Also
        --------
        Dataset.xindexes

        """
        return self.xindexes.to_pandas_indexes()

    @property
    def xindexes(self) -> Indexes[Index]:
        """Mapping of :py:class:`~xarray.indexes.Index` objects
        used for label based indexing.
        """
        return Indexes(self._indexes, {k: self._variables[k] for k in self._indexes})

    @property
    def coords(self) -> DatasetCoordinates:
        """Mapping of :py:class:`~xarray.DataArray` objects corresponding to
        coordinate variables.

        See Also
        --------
        Coordinates
        """
        return DatasetCoordinates(self)

    @property
    def data_vars(self) -> DataVariables:
        """Dictionary of DataArray objects corresponding to data variables"""
        return DataVariables(self)

    def set_coords(self, names: Hashable | Iterable[Hashable]) -> Self:
        """Given names of one or more variables, set them as coordinates

        Parameters
        ----------
        names : hashable or iterable of hashable
            Name(s) of variables in this dataset to convert into coordinates.

        Examples
        --------
        >>> dataset = xr.Dataset(
        ...     {
        ...         "pressure": ("time", [1.013, 1.2, 3.5]),
        ...         "time": pd.date_range("2023-01-01", periods=3),
        ...     }
        ... )
        >>> dataset
        <xarray.Dataset> Size: 48B
        Dimensions:   (time: 3)
        Coordinates:
          * time      (time) datetime64[ns] 24B 2023-01-01 2023-01-02 2023-01-03
        Data variables:
            pressure  (time) float64 24B 1.013 1.2 3.5

        >>> dataset.set_coords("pressure")
        <xarray.Dataset> Size: 48B
        Dimensions:   (time: 3)
        Coordinates:
            pressure  (time) float64 24B 1.013 1.2 3.5
          * time      (time) datetime64[ns] 24B 2023-01-01 2023-01-02 2023-01-03
        Data variables:
            *empty*

        On calling ``set_coords`` , these data variables are converted to coordinates, as shown in the final dataset.

        Returns
        -------
        Dataset

        See Also
        --------
        Dataset.swap_dims
        Dataset.assign_coords
        """
        # TODO: allow inserting new coordinates with this method, like
        # DataFrame.set_index?
        # nb. check in self._variables, not self.data_vars to insure that the
        # operation is idempotent
        if isinstance(names, str) or not isinstance(names, Iterable):
            names = [names]
        else:
            names = list(names)
        self._assert_all_in_dataset(names)
        obj = self.copy()
        obj._coord_names.update(names)
        return obj

    def reset_coords(
        self,
        names: Dims = None,
        drop: bool = False,
    ) -> Self:
        """Given names of coordinates, reset them to become variables

        Parameters
        ----------
        names : str, Iterable of Hashable or None, optional
            Name(s) of non-index coordinates in this dataset to reset into
            variables. By default, all non-index coordinates are reset.
        drop : bool, default: False
            If True, remove coordinates instead of converting them into
            variables.

        Examples
        --------
        >>> dataset = xr.Dataset(
        ...     {
        ...         "temperature": (
        ...             ["time", "lat", "lon"],
        ...             [[[25, 26], [27, 28]], [[29, 30], [31, 32]]],
        ...         ),
        ...         "precipitation": (
        ...             ["time", "lat", "lon"],
        ...             [[[0.5, 0.8], [0.2, 0.4]], [[0.3, 0.6], [0.7, 0.9]]],
        ...         ),
        ...     },
        ...     coords={
        ...         "time": pd.date_range(start="2023-01-01", periods=2),
        ...         "lat": [40, 41],
        ...         "lon": [-80, -79],
        ...         "altitude": 1000,
        ...     },
        ... )

        # Dataset before resetting coordinates

        >>> dataset
        <xarray.Dataset> Size: 184B
        Dimensions:        (time: 2, lat: 2, lon: 2)
        Coordinates:
          * time           (time) datetime64[ns] 16B 2023-01-01 2023-01-02
          * lat            (lat) int64 16B 40 41
          * lon            (lon) int64 16B -80 -79
            altitude       int64 8B 1000
        Data variables:
            temperature    (time, lat, lon) int64 64B 25 26 27 28 29 30 31 32
            precipitation  (time, lat, lon) float64 64B 0.5 0.8 0.2 0.4 0.3 0.6 0.7 0.9

        # Reset the 'altitude' coordinate

        >>> dataset_reset = dataset.reset_coords("altitude")

        # Dataset after resetting coordinates

        >>> dataset_reset
        <xarray.Dataset> Size: 184B
        Dimensions:        (time: 2, lat: 2, lon: 2)
        Coordinates:
          * time           (time) datetime64[ns] 16B 2023-01-01 2023-01-02
          * lat            (lat) int64 16B 40 41
          * lon            (lon) int64 16B -80 -79
        Data variables:
            temperature    (time, lat, lon) int64 64B 25 26 27 28 29 30 31 32
            precipitation  (time, lat, lon) float64 64B 0.5 0.8 0.2 0.4 0.3 0.6 0.7 0.9
            altitude       int64 8B 1000

        Returns
        -------
        Dataset

        See Also
        --------
        Dataset.set_coords
        """
        if names is None:
            names = self._coord_names - set(self._indexes)
        else:
            if isinstance(names, str) or not isinstance(names, Iterable):
                names = [names]
            else:
                names = list(names)
            self._assert_all_in_dataset(names)
            bad_coords = set(names) & set(self._indexes)
            if bad_coords:
                raise ValueError(
                    f"cannot remove index coordinates with reset_coords: {bad_coords}"
                )
        obj = self.copy()
        obj._coord_names.difference_update(names)
        if drop:
            for name in names:
                del obj._variables[name]
        return obj

    def dump_to_store(self, store: AbstractDataStore, **kwargs) -> None:
        """Store dataset contents to a backends.*DataStore object."""
        from xarray.backends.api import dump_to_store

        # TODO: rename and/or cleanup this method to make it more consistent
        # with to_netcdf()
        dump_to_store(self, store, **kwargs)

    # path=None writes to bytes
    @overload
    def to_netcdf(
        self,
        path: None = None,
        mode: NetcdfWriteModes = "w",
        format: T_NetcdfTypes | None = None,
        group: str | None = None,
        engine: T_NetcdfEngine | None = None,
        encoding: Mapping[Any, Mapping[str, Any]] | None = None,
        unlimited_dims: Iterable[Hashable] | None = None,
        compute: bool = True,
        invalid_netcdf: bool = False,
        auto_complex: bool | None = None,
    ) -> memoryview: ...

    # compute=False returns dask.Delayed
    @overload
    def to_netcdf(
        self,
        path: str | PathLike,
        mode: NetcdfWriteModes = "w",
        format: T_NetcdfTypes | None = None,
        group: str | None = None,
        engine: T_NetcdfEngine | None = None,
        encoding: Mapping[Any, Mapping[str, Any]] | None = None,
        unlimited_dims: Iterable[Hashable] | None = None,
        *,
        compute: Literal[False],
        invalid_netcdf: bool = False,
        auto_complex: bool | None = None,
    ) -> Delayed: ...

    # default return None
    @overload
    def to_netcdf(
        self,
        path: str | PathLike | io.IOBase,
        mode: NetcdfWriteModes = "w",
        format: T_NetcdfTypes | None = None,
        group: str | None = None,
        engine: T_NetcdfEngine | None = None,
        encoding: Mapping[Any, Mapping[str, Any]] | None = None,
        unlimited_dims: Iterable[Hashable] | None = None,
        compute: Literal[True] = True,
        invalid_netcdf: bool = False,
        auto_complex: bool | None = None,
    ) -> None: ...

    # if compute cannot be evaluated at type check time
    # we may get back either Delayed or None
    @overload
    def to_netcdf(
        self,
        path: str | PathLike,
        mode: NetcdfWriteModes = "w",
        format: T_NetcdfTypes | None = None,
        group: str | None = None,
        engine: T_NetcdfEngine | None = None,
        encoding: Mapping[Any, Mapping[str, Any]] | None = None,
        unlimited_dims: Iterable[Hashable] | None = None,
        compute: bool = True,
        invalid_netcdf: bool = False,
        auto_complex: bool | None = None,
    ) -> Delayed | None: ...

    def to_netcdf(
        self,
        path: str | PathLike | io.IOBase | None = None,
        mode: NetcdfWriteModes = "w",
        format: T_NetcdfTypes | None = None,
        group: str | None = None,
        engine: T_NetcdfEngine | None = None,
        encoding: Mapping[Any, Mapping[str, Any]] | None = None,
        unlimited_dims: Iterable[Hashable] | None = None,
        compute: bool = True,
        invalid_netcdf: bool = False,
        auto_complex: bool | None = None,
    ) -> memoryview | Delayed | None:
        """Write dataset contents to a netCDF file.

        Parameters
        ----------
        path : str, path-like, file-like or None, optional
            Path to which to save this datatree, or a file-like object to write
            it to (which must support read and write and be seekable) or None
            (default) to return in-memory bytes as a memoryview.
        mode : {"w", "a"}, default: "w"
            Write ('w') or append ('a') mode. If mode='w', any existing file at
            this location will be overwritten. If mode='a', existing variables
            will be overwritten.
        format : {"NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_64BIT", \
                  "NETCDF3_CLASSIC"}, optional
            File format for the resulting netCDF file:

            * NETCDF4: Data is stored in an HDF5 file, using netCDF4 API
              features.
            * NETCDF4_CLASSIC: Data is stored in an HDF5 file, using only
              netCDF 3 compatible API features.
            * NETCDF3_64BIT: 64-bit offset version of the netCDF 3 file format,
              which fully supports 2+ GB files, but is only compatible with
              clients linked against netCDF version 3.6.0 or later.
            * NETCDF3_CLASSIC: The classic netCDF 3 file format. It does not
              handle 2+ GB files very well.

            All formats are supported by the netCDF4-python library.
            scipy.io.netcdf only supports the last two formats.

            The default format is NETCDF4 if you are saving a file to disk and
            have the netCDF4-python library available. Otherwise, xarray falls
            back to using scipy to write netCDF files and defaults to the
            NETCDF3_64BIT format (scipy does not support netCDF4).
        group : str, optional
            Path to the netCDF4 group in the given file to open (only works for
            format='NETCDF4'). The group(s) will be created if necessary.
        engine : {"netcdf4", "scipy", "h5netcdf"}, optional
            Engine to use when writing netCDF files. If not provided, the
            default engine is chosen based on available dependencies, with a
            preference for 'netcdf4' if writing to a file on disk.
        encoding : dict, optional
            Nested dictionary with variable names as keys and dictionaries of
            variable specific encodings as values, e.g.,
            ``{"my_variable": {"dtype": "int16", "scale_factor": 0.1,
            "zlib": True}, ...}``.
            If ``encoding`` is specified the original encoding of the variables of
            the dataset is ignored.

            The `h5netcdf` engine supports both the NetCDF4-style compression
            encoding parameters ``{"zlib": True, "complevel": 9}`` and the h5py
            ones ``{"compression": "gzip", "compression_opts": 9}``.
            This allows using any compression plugin installed in the HDF5
            library, e.g. LZF.

        unlimited_dims : iterable of hashable, optional
            Dimension(s) that should be serialized as unlimited dimensions.
            By default, no dimensions are treated as unlimited dimensions.
            Note that unlimited_dims may also be set via
            ``dataset.encoding["unlimited_dims"]``.
        compute: bool, default: True
            If true compute immediately, otherwise return a
            ``dask.delayed.Delayed`` object that can be computed later.
        invalid_netcdf: bool, default: False
            Only valid along with ``engine="h5netcdf"``. If True, allow writing
            hdf5 files which are invalid netcdf as described in
            https://github.com/h5netcdf/h5netcdf.

        Returns
        -------
            * ``memoryview`` if path is None
            * ``dask.delayed.Delayed`` if compute is False
            * ``None`` otherwise

        See Also
        --------
        DataArray.to_netcdf
        """
        if encoding is None:
            encoding = {}
        from xarray.backends.api import to_netcdf

        return to_netcdf(  # type: ignore[return-value]  # mypy cannot resolve the overloads:(
            self,
            path,
            mode=mode,
            format=format,
            group=group,
            engine=engine,
            encoding=encoding,
            unlimited_dims=unlimited_dims,
            compute=compute,
            multifile=False,
            invalid_netcdf=invalid_netcdf,
            auto_complex=auto_complex,
        )

    # compute=True (default) returns ZarrStore
    @overload
    def to_zarr(
        self,
        store: MutableMapping | str | PathLike[str] | None = None,
        chunk_store: MutableMapping | str | PathLike | None = None,
        mode: ZarrWriteModes | None = None,
        synchronizer=None,
        group: str | None = None,
        encoding: Mapping | None = None,
        *,
        compute: Literal[True] = True,
        consolidated: bool | None = None,
        append_dim: Hashable | None = None,
        region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None = None,
        safe_chunks: bool = True,
        align_chunks: bool = False,
        storage_options: dict[str, str] | None = None,
        zarr_version: int | None = None,
        zarr_format: int | None = None,
        write_empty_chunks: bool | None = None,
        chunkmanager_store_kwargs: dict[str, Any] | None = None,
    ) -> ZarrStore: ...

    # compute=False returns dask.Delayed
    @overload
    def to_zarr(
        self,
        store: MutableMapping | str | PathLike[str] | None = None,
        chunk_store: MutableMapping | str | PathLike | None = None,
        mode: ZarrWriteModes | None = None,
        synchronizer=None,
        group: str | None = None,
        encoding: Mapping | None = None,
        *,
        compute: Literal[False],
        consolidated: bool | None = None,
        append_dim: Hashable | None = None,
        region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None = None,
        safe_chunks: bool = True,
        align_chunks: bool = False,
        storage_options: dict[str, str] | None = None,
        zarr_version: int | None = None,
        zarr_format: int | None = None,
        write_empty_chunks: bool | None = None,
        chunkmanager_store_kwargs: dict[str, Any] | None = None,
    ) -> Delayed: ...

    def to_zarr(
        self,
        store: MutableMapping | str | PathLike[str] | None = None,
        chunk_store: MutableMapping | str | PathLike | None = None,
        mode: ZarrWriteModes | None = None,
        synchronizer=None,
        group: str | None = None,
        encoding: Mapping | None = None,
        *,
        compute: bool = True,
        consolidated: bool | None = None,
        append_dim: Hashable | None = None,
        region: Mapping[str, slice | Literal["auto"]] | Literal["auto"] | None = None,
        safe_chunks: bool = True,
        align_chunks: bool = False,
        storage_options: dict[str, str] | None = None,
        zarr_version: int | None = None,
        zarr_format: int | None = None,
        write_empty_chunks: bool | None = None,
        chunkmanager_store_kwargs: dict[str, Any] | None = None,
    ) -> ZarrStore | Delayed:
        """Write dataset contents to a zarr group.

        Zarr chunks are determined in the following way:

        - From the ``chunks`` attribute in each variable's ``encoding``
          (can be set via `Dataset.chunk`).
        - If the variable is a Dask array, from the dask chunks
        - If neither Dask chunks nor encoding chunks are present, chunks will
          be determined automatically by Zarr
        - If both Dask chunks and encoding chunks are present, encoding chunks
          will be used, provided that there is a many-to-one relationship between
          encoding chunks and dask chunks (i.e. Dask chunks are bigger than and
          evenly divide encoding chunks); otherwise raise a ``ValueError``.
          This restriction ensures that no synchronization / locks are required
          when writing. To disable this restriction, use ``safe_chunks=False``.

        Parameters
        ----------
        store : MutableMapping, str or path-like, optional
            Store or path to directory in local or remote file system.
        chunk_store : MutableMapping, str or path-like, optional
            Store or path to directory in local or remote file system only for Zarr
            array chunks. Requires zarr-python v2.4.0 or later.
        mode : {"w", "w-", "a", "a-", r+", None}, optional
            Persistence mode: "w" means create (overwrite if exists);
            "w-" means create (fail if exists);
            "a" means override all existing variables including dimension coordinates (create if does not exist);
            "a-" means only append those variables that have ``append_dim``.
            "r+" means modify existing array *values* only (raise an error if
            any metadata or shapes would change).
            The default mode is "a" if ``append_dim`` is set. Otherwise, it is
            "r+" if ``region`` is set and ``w-`` otherwise.
        synchronizer : object, optional
            Zarr array synchronizer.
        group : str, optional
            Group path. (a.k.a. `path` in zarr terminology.)
        encoding : dict, optional
            Nested dictionary with variable names as keys and dictionaries of
            variable specific encodings as values, e.g.,
            ``{"my_variable": {"dtype": "int16", "scale_factor": 0.1,}, ...}``
        compute : bool, default: True
            If True write array data immediately, otherwise return a
            ``dask.delayed.Delayed`` object that can be computed to write
            array data later. Metadata is always updated eagerly.
        consolidated : bool, optional
            If True, apply :func:`zarr.convenience.consolidate_metadata`
            after writing metadata and read existing stores with consolidated
            metadata; if False, do not. The default (`consolidated=None`) means
            write consolidated metadata and attempt to read consolidated
            metadata for existing stores (falling back to non-consolidated).

            When the experimental ``zarr_version=3``, ``consolidated`` must be
            either be ``None`` or ``False``.
        append_dim : hashable, optional
            If set, the dimension along which the data will be appended. All
            other dimensions on overridden variables must remain the same size.
        region : dict or "auto", optional
            Optional mapping from dimension names to either a) ``"auto"``, or b) integer
            slices, indicating the region of existing zarr array(s) in which to write
            this dataset's data.

            If ``"auto"`` is provided the existing store will be opened and the region
            inferred by matching indexes. ``"auto"`` can be used as a single string,
            which will automatically infer the region for all dimensions, or as
            dictionary values for specific dimensions mixed together with explicit
            slices for other dimensions.

            Alternatively integer slices can be provided; for example, ``{'x': slice(0,
            1000), 'y': slice(10000, 11000)}`` would indicate that values should be
            written to the region ``0:1000`` along ``x`` and ``10000:11000`` along
            ``y``.

            Two restrictions apply to the use of ``region``:

            - If ``region`` is set, _all_ variables in a dataset must have at
              least one dimension in common with the region. Other variables
              should be written in a separate single call to ``to_zarr()``.
            - Dimensions cannot be included in both ``region`` and
              ``append_dim`` at the same time. To create empty arrays to fill
              in with ``region``, use a separate call to ``to_zarr()`` with
              ``compute=False``. See "Modifying existing Zarr stores" in
              the reference documentation for full details.

            Users are expected to ensure that the specified region aligns with
            Zarr chunk boundaries, and that dask chunks are also aligned.
            Xarray makes limited checks that these multiple chunk boundaries line up.
            It is possible to write incomplete chunks and corrupt the data with this
            option if you are not careful.
        safe_chunks : bool, default: True
            If True, only allow writes to when there is a many-to-one relationship
            between Zarr chunks (specified in encoding) and Dask chunks.
            Set False to override this restriction; however, data may become corrupted
            if Zarr arrays are written in parallel. This option may be useful in combination
            with ``compute=False`` to initialize a Zarr from an existing
            Dataset with arbitrary chunk structure.
            In addition to the many-to-one relationship validation, it also detects partial
            chunks writes when using the region parameter,
            these partial chunks are considered unsafe in the mode "r+" but safe in
            the mode "a".
            Note: Even with these validations it can still be unsafe to write
            two or more chunked arrays in the same location in parallel if they are
            not writing in independent regions, for those cases it is better to use
            a synchronizer.
        align_chunks: bool, default False
            If True, rechunks the Dask array to align with Zarr chunks before writing.
            This ensures each Dask chunk maps to one or more contiguous Zarr chunks,
            which avoids race conditions.
            Internally, the process sets safe_chunks=False and tries to preserve
            the original Dask chunking as much as possible.
            Note: While this alignment avoids write conflicts stemming from chunk
            boundary misalignment, it does not protect against race conditions
            if multiple uncoordinated processes write to the same
            Zarr array concurrently.
        storage_options : dict, optional
            Any additional parameters for the storage backend (ignored for local
            paths).
        zarr_version : int or None, optional

            .. deprecated:: 2024.9.1
            Use ``zarr_format`` instead.

        zarr_format : int or None, optional
            The desired zarr format to target (currently 2 or 3). The default
            of None will attempt to determine the zarr version from ``store`` when
            possible, otherwise defaulting to the default version used by
            the zarr-python library installed.
        write_empty_chunks : bool or None, optional
            If True, all chunks will be stored regardless of their
            contents. If False, each chunk is compared to the array's fill value
            prior to storing. If a chunk is uniformly equal to the fill value, then
            that chunk is not be stored, and the store entry for that chunk's key
            is deleted. This setting enables sparser storage, as only chunks with
            non-fill-value data are stored, at the expense of overhead associated
            with checking the data of each chunk. If None (default) fall back to
            specification(s) in ``encoding`` or Zarr defaults. A ``ValueError``
            will be raised if the value of this (if not None) differs with
            ``encoding``.
        chunkmanager_store_kwargs : dict, optional
            Additional keyword arguments passed on to the `ChunkManager.store` method used to store
            chunked arrays. For example for a dask array additional kwargs will be passed eventually to
            :py:func:`dask.array.store()`. Experimental API that should not be relied upon.

        Returns
        -------
            * ``dask.delayed.Delayed`` if compute is False
            * ZarrStore otherwise

        References
        ----------
        https://zarr.readthedocs.io/

        Notes
        -----
        Zarr chunking behavior:
            If chunks are found in the encoding argument or attribute
            corresponding to any DataArray, those chunks are used.
            If a DataArray is a dask array, it is written with those chunks.
            If not other chunks are found, Zarr uses its own heuristics to
            choose automatic chunk sizes.

        encoding:
            The encoding attribute (if exists) of the DataArray(s) will be
            used. Override any existing encodings by providing the ``encoding`` kwarg.

        ``fill_value`` handling:
            There exists a subtlety in interpreting zarr's ``fill_value`` property. For zarr v2 format
            arrays, ``fill_value`` is *always* interpreted as an invalid value similar to the ``_FillValue`` attribute
            in CF/netCDF. For Zarr v3 format arrays, only an explicit ``_FillValue`` attribute will be used
            to mask the data if requested using ``mask_and_scale=True``. See this `Github issue <https://github.com/pydata/xarray/issues/5475>`_
            for more.

        See Also
        --------
        :ref:`io.zarr`
            The I/O user guide, with more details and examples.
        """
        from xarray.backends.api import to_zarr

        return to_zarr(  # type: ignore[call-overload,misc]
            self,
            store=store,
            chunk_store=chunk_store,
            storage_options=storage_options,
            mode=mode,
            synchronizer=synchronizer,
            group=group,
            encoding=encoding,
            compute=compute,
            consolidated=consolidated,
            append_dim=append_dim,
            region=region,
            safe_chunks=safe_chunks,
            zarr_version=zarr_version,
            zarr_format=zarr_format,
            write_empty_chunks=write_empty_chunks,
            chunkmanager_store_kwargs=chunkmanager_store_kwargs,
        )

    def __repr__(self) -> str:
        return formatting.dataset_repr(self)

    def _repr_html_(self) -> str:
        if OPTIONS["display_style"] == "text":
            return f"<pre>{escape(repr(self))}</pre>"
        return formatting_html.dataset_repr(self)

    def info(self, buf: IO | None = None) -> None:
        """
        Concise summary of a Dataset variables and attributes.

        Parameters
        ----------
        buf : file-like, default: sys.stdout
            writable buffer

        See Also
        --------
        pandas.DataFrame.assign
        ncdump : netCDF's ncdump
        """
        if buf is None:  # pragma: no cover
            buf = sys.stdout

        lines = [
            "xarray.Dataset {",
            "dimensions:",
        ]
        for name, size in self.sizes.items():
            lines.append(f"\t{name} = {size} ;")
        lines.append("\nvariables:")
        for name, da in self.variables.items():
            dims = ", ".join(map(str, da.dims))
            lines.append(f"\t{da.dtype} {name}({dims}) ;")
            for k, v in da.attrs.items():
                lines.append(f"\t\t{name}:{k} = {v} ;")
        lines.append("\n// global attributes:")
        for k, v in self.attrs.items():
            lines.append(f"\t:{k} = {v} ;")
        lines.append("}")

        buf.write("\n".join(lines))

    @property
    def chunks(self) -> Mapping[Hashable, tuple[int, ...]]:
        """
        Mapping from dimension names to block lengths for this dataset's data.

        If this dataset does not contain chunked arrays, the mapping will be empty.

        Cannot be modified directly, but can be modified by calling .chunk().

        Same as Dataset.chunksizes, but maintained for backwards compatibility.

        See Also
        --------
        Dataset.chunk
        Dataset.chunksizes
        xarray.unify_chunks
        """
        return get_chunksizes(self.variables.values())

    @property
    def chunksizes(self) -> Mapping[Hashable, tuple[int, ...]]:
        """
        Mapping from dimension names to block lengths for this dataset's data.

        If this dataset does not contain chunked arrays, the mapping will be empty.

        Cannot be modified directly, but can be modified by calling .chunk().

        Same as Dataset.chunks.

        See Also
        --------
        Dataset.chunk
        Dataset.chunks
        xarray.unify_chunks
        """
        return get_chunksizes(self.variables.values())

    def chunk(
        self,
        chunks: T_ChunksFreq = {},  # noqa: B006  # {} even though it's technically unsafe, is being used intentionally here (#4667)
        name_prefix: str = "xarray-",
        token: str | None = None,
        lock: bool = False,
        inline_array: bool = False,
        chunked_array_type: str | ChunkManagerEntrypoint | None = None,
        from_array_kwargs=None,
        **chunks_kwargs: T_ChunkDimFreq,
    ) -> Self:
        """Coerce all arrays in this dataset into dask arrays with the given
        chunks.

        Non-dask arrays in this dataset will be converted to dask arrays. Dask
        arrays will be rechunked to the given chunk sizes.

        If neither chunks is not provided for one or more dimensions, chunk
        sizes along that dimension will not be updated; non-dask arrays will be
        converted into dask arrays with a single block.

        Along datetime-like dimensions, a :py:class:`Resampler` object
        (e.g. :py:class:`groupers.TimeResampler` or :py:class:`groupers.SeasonResampler`)
        is also accepted.

        Parameters
        ----------
        chunks : int, tuple of int, "auto" or mapping of hashable to int or a Resampler, optional
            Chunk sizes along each dimension, e.g., ``5``, ``"auto"``, or
            ``{"x": 5, "y": 5}`` or ``{"x": 5, "time": TimeResampler(freq="YE")}`` or
            ``{"time": SeasonResampler(["DJF", "MAM", "JJA", "SON"])}``.
        name_prefix : str, default: "xarray-"
            Prefix for the name of any new dask arrays.
        token : str, optional
            Token uniquely identifying this dataset.
        lock : bool, default: False
            Passed on to :py:func:`dask.array.from_array`, if the array is not
            already as dask array.
        inline_array: bool, default: False
            Passed on to :py:func:`dask.array.from_array`, if the array is not
            already as dask array.
        chunked_array_type: str, optional
            Which chunked array type to coerce this datasets' arrays to.
            Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEntryPoint` system.
            Experimental API that should not be relied upon.
        from_array_kwargs: dict, optional
            Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
            chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
            For example, with dask as the default chunked array type, this method would pass additional kwargs
            to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
        **chunks_kwargs : {dim: chunks, ...}, optional
            The keyword arguments form of ``chunks``.
            One of chunks or chunks_kwargs must be provided

        Returns
        -------
        chunked : xarray.Dataset

        See Also
        --------
        Dataset.chunks
        Dataset.chunksizes
        xarray.unify_chunks
        dask.array.from_array
        """
        from xarray.groupers import Resampler

        if chunks is None and not chunks_kwargs:
            warnings.warn(
                "None value for 'chunks' is deprecated. "
                "It will raise an error in the future. Use instead '{}'",
                category=DeprecationWarning,
                stacklevel=2,
            )
            chunks = {}
        chunks_mapping: Mapping[Any, Any]
        if not isinstance(chunks, Mapping) and chunks is not None:
            if isinstance(chunks, tuple | list):
                utils.emit_user_level_warning(
                    "Supplying chunks as dimension-order tuples is deprecated. "
                    "It will raise an error in the future. Instead use a dict with dimensions as keys.",
                    category=DeprecationWarning,
                )
            chunks_mapping = dict.fromkeys(self.dims, chunks)
        else:
            chunks_mapping = either_dict_or_kwargs(chunks, chunks_kwargs, "chunk")

        bad_dims = chunks_mapping.keys() - self.sizes.keys()
        if bad_dims:
            raise ValueError(
                f"chunks keys {tuple(bad_dims)} not found in data dimensions {tuple(self.sizes.keys())}"
            )

        def _resolve_resampler(name: Hashable, resampler: Resampler) -> tuple[int, ...]:
            variable = self._variables.get(name, None)
            if variable is None:
                raise ValueError(
                    f"Cannot chunk by resampler {resampler!r} for virtual variable {name!r}."
                )
            if variable.ndim != 1:
                raise ValueError(
                    f"chunks={resampler!r} only supported for 1D variables. "
                    f"Received variable {name!r} with {variable.ndim} dimensions instead."
                )
            newchunks = resampler.compute_chunks(variable, dim=name)
            if sum(newchunks) != variable.shape[0]:
                raise ValueError(
                    f"Logic bug in rechunking variable {name!r} using {resampler!r}. "
                    "New chunks tuple does not match size of data. Please open an issue."
                )
            return newchunks

        chunks_mapping_ints: Mapping[Any, T_ChunkDim] = {
            name: (
                _resolve_resampler(name, chunks)
                if isinstance(chunks, Resampler)
                else chunks
            )
            for name, chunks in chunks_mapping.items()
        }

        chunkmanager = guess_chunkmanager(chunked_array_type)
        if from_array_kwargs is None:
            from_array_kwargs = {}

        variables = {
            k: _maybe_chunk(
                k,
                v,
                chunks_mapping_ints,
                token,
                lock,
                name_prefix,
                inline_array=inline_array,
                chunked_array_type=chunkmanager,
                from_array_kwargs=from_array_kwargs.copy(),
            )
            for k, v in self.variables.items()
        }
        return self._replace(variables)

    def _validate_indexers(
        self, indexers: Mapping[Any, Any], missing_dims: ErrorOptionsWithWarn = "raise"
    ) -> Iterator[tuple[Hashable, int | slice | np.ndarray | Variable]]:
        """Here we make sure
        + indexer has a valid keys
        + indexer is in a valid data type
        + string indexers are cast to the appropriate date type if the
          associated index is a DatetimeIndex or CFTimeIndex
        """
        from xarray.core.dataarray import DataArray

        indexers = drop_dims_from_indexers(indexers, self.dims, missing_dims)

        # all indexers should be int, slice, np.ndarrays, or Variable
        for k, v in indexers.items():
            if isinstance(v, int | slice | Variable) and not isinstance(v, bool):
                yield k, v
            elif isinstance(v, DataArray):
                yield k, v.variable
            elif isinstance(v, tuple):
                yield k, as_variable(v)
            elif isinstance(v, Dataset):
                raise TypeError("cannot use a Dataset as an indexer")
            elif isinstance(v, Sequence) and len(v) == 0:
                yield k, np.empty((0,), dtype="int64")
            else:
                if not is_duck_array(v):
                    v = np.asarray(v)

                if v.dtype.kind in "US":
                    index = self._indexes[k].to_pandas_index()
                    if isinstance(index, pd.DatetimeIndex):
                        v = duck_array_ops.astype(v, dtype="datetime64[ns]")
                    elif isinstance(index, CFTimeIndex):
                        v = _parse_array_of_cftime_strings(v, index.date_type)

                if v.ndim > 1:
                    raise IndexError(
                        "Unlabeled multi-dimensional array cannot be "
                        f"used for indexing: {k}"
                    )
                yield k, v

    def _validate_interp_indexers(
        self, indexers: Mapping[Any, Any]
    ) -> Iterator[tuple[Hashable, Variable]]:
        """Variant of _validate_indexers to be used for interpolation"""
        for k, v in self._validate_indexers(indexers):
            if isinstance(v, Variable):
                yield k, v
            elif is_scalar(v):
                yield k, Variable((), v, attrs=self.coords[k].attrs)
            elif isinstance(v, np.ndarray):
                yield k, Variable(dims=(k,), data=v, attrs=self.coords[k].attrs)
            else:
                raise TypeError(type(v))

    def _get_indexers_coords_and_indexes(self, indexers):
        """Extract coordinates and indexes from indexers.

        Only coordinate with a name different from any of self.variables will
        be attached.
        """
        from xarray.core.dataarray import DataArray

        coords_list = []
        for k, v in indexers.items():
            if isinstance(v, DataArray):
                if v.dtype.kind == "b":
                    if v.ndim != 1:  # we only support 1-d boolean array
                        raise ValueError(
                            f"{v.ndim:d}d-boolean array is used for indexing along "
                            f"dimension {k!r}, but only 1d boolean arrays are "
                            "supported."
                        )
                    # Make sure in case of boolean DataArray, its
                    # coordinate also should be indexed.
                    v_coords = v[v.values.nonzero()[0]].coords
                else:
                    v_coords = v.coords
                coords_list.append(v_coords)

        # we don't need to call align() explicitly or check indexes for
        # alignment, because merge_variables already checks for exact alignment
        # between dimension coordinates
        coords, indexes = merge_coordinates_without_align(coords_list)
        assert_coordinate_consistent(self, coords)

        # silently drop the conflicted variables.
        attached_coords = {k: v for k, v in coords.items() if k not in self._variables}
        attached_indexes = {
            k: v for k, v in indexes.items() if k not in self._variables
        }
        return attached_coords, attached_indexes

    def isel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        drop: bool = False,
        missing_dims: ErrorOptionsWithWarn = "raise",
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new dataset with each array indexed along the specified
        dimension(s).

        This method selects values from each array using its `__getitem__`
        method, except this method does not require knowing the order of
        each array's dimensions.

        Parameters
        ----------
        indexers : dict, optional
            A dict with keys matching dimensions and values given
            by integers, slice objects or arrays.
            indexer can be a integer, slice, array-like or DataArray.
            If DataArrays are passed as indexers, xarray-style indexing will be
            carried out. See :ref:`indexing` for the details.
            One of indexers or indexers_kwargs must be provided.
        drop : bool, default: False
            If ``drop=True``, drop coordinates variables indexed by integers
            instead of making them scalar.
        missing_dims : {"raise", "warn", "ignore"}, default: "raise"
            What to do if dimensions that should be selected from are not present in the
            Dataset:
            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        **indexers_kwargs : {dim: indexer, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Returns
        -------
        obj : Dataset
            A new Dataset with the same contents as this dataset, except each
            array and dimension is indexed by the appropriate indexers.
            If indexer DataArrays have coordinates that do not conflict with
            this object, then these coordinates will be attached.
            In general, each array's data will be a view of the array's data
            in this dataset, unless vectorized indexing was triggered by using
            an array indexer, in which case the data will be a copy.

        Examples
        --------

        >>> dataset = xr.Dataset(
        ...     {
        ...         "math_scores": (
        ...             ["student", "test"],
        ...             [[90, 85, 92], [78, 80, 85], [95, 92, 98]],
        ...         ),
        ...         "english_scores": (
        ...             ["student", "test"],
        ...             [[88, 90, 92], [75, 82, 79], [93, 96, 91]],
        ...         ),
        ...     },
        ...     coords={
        ...         "student": ["Alice", "Bob", "Charlie"],
        ...         "test": ["Test 1", "Test 2", "Test 3"],
        ...     },
        ... )

        # A specific element from the dataset is selected

        >>> dataset.isel(student=1, test=0)
        <xarray.Dataset> Size: 68B
        Dimensions:         ()
        Coordinates:
            student         <U7 28B 'Bob'
            test            <U6 24B 'Test 1'
        Data variables:
            math_scores     int64 8B 78
            english_scores  int64 8B 75

        # Indexing with a slice using isel

        >>> slice_of_data = dataset.isel(student=slice(0, 2), test=slice(0, 2))
        >>> slice_of_data
        <xarray.Dataset> Size: 168B
        Dimensions:         (student: 2, test: 2)
        Coordinates:
          * student         (student) <U7 56B 'Alice' 'Bob'
          * test            (test) <U6 48B 'Test 1' 'Test 2'
        Data variables:
            math_scores     (student, test) int64 32B 90 85 78 80
            english_scores  (student, test) int64 32B 88 90 75 82

        >>> index_array = xr.DataArray([0, 2], dims="student")
        >>> indexed_data = dataset.isel(student=index_array)
        >>> indexed_data
        <xarray.Dataset> Size: 224B
        Dimensions:         (student: 2, test: 3)
        Coordinates:
          * student         (student) <U7 56B 'Alice' 'Charlie'
          * test            (test) <U6 72B 'Test 1' 'Test 2' 'Test 3'
        Data variables:
            math_scores     (student, test) int64 48B 90 85 92 95 92 98
            english_scores  (student, test) int64 48B 88 90 92 93 96 91

        See Also
        --------
        :func:`Dataset.sel <Dataset.sel>`
        :func:`DataArray.isel <DataArray.isel>`

        :doc:`xarray-tutorial:intermediate/indexing/indexing`
            Tutorial material on indexing with Xarray objects

        :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
            Tutorial material on basics of indexing

        """
        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "isel")
        if any(is_fancy_indexer(idx) for idx in indexers.values()):
            return self._isel_fancy(indexers, drop=drop, missing_dims=missing_dims)

        # Much faster algorithm for when all indexers are ints, slices, one-dimensional
        # lists, or zero or one-dimensional np.ndarray's
        indexers = drop_dims_from_indexers(indexers, self.dims, missing_dims)

        variables = {}
        dims: dict[Hashable, int] = {}
        coord_names = self._coord_names.copy()

        indexes, index_variables = isel_indexes(self.xindexes, indexers)

        for name, var in self._variables.items():
            # preserve variable order
            if name in index_variables:
                var = index_variables[name]
            else:
                var_indexers = {k: v for k, v in indexers.items() if k in var.dims}
                if var_indexers:
                    var = var.isel(var_indexers)
                    if drop and var.ndim == 0 and name in coord_names:
                        coord_names.remove(name)
                        continue
            variables[name] = var
            dims.update(zip(var.dims, var.shape, strict=True))

        return self._construct_direct(
            variables=variables,
            coord_names=coord_names,
            dims=dims,
            attrs=self._attrs,
            indexes=indexes,
            encoding=self._encoding,
            close=self._close,
        )

    def _isel_fancy(
        self,
        indexers: Mapping[Any, Any],
        *,
        drop: bool,
        missing_dims: ErrorOptionsWithWarn = "raise",
    ) -> Self:
        valid_indexers = dict(self._validate_indexers(indexers, missing_dims))

        variables: dict[Hashable, Variable] = {}
        indexes, index_variables = isel_indexes(self.xindexes, valid_indexers)

        for name, var in self.variables.items():
            if name in index_variables:
                new_var = index_variables[name]
            else:
                var_indexers = {
                    k: v for k, v in valid_indexers.items() if k in var.dims
                }
                if var_indexers:
                    new_var = var.isel(indexers=var_indexers)
                    # drop scalar coordinates
                    # https://github.com/pydata/xarray/issues/6554
                    if name in self.coords and drop and new_var.ndim == 0:
                        continue
                else:
                    new_var = var.copy(deep=False)
                if name not in indexes:
                    new_var = new_var.to_base_variable()
            variables[name] = new_var

        coord_names = self._coord_names & variables.keys()
        selected = self._replace_with_new_dims(variables, coord_names, indexes)

        # Extract coordinates from indexers
        coord_vars, new_indexes = selected._get_indexers_coords_and_indexes(indexers)
        variables.update(coord_vars)
        indexes.update(new_indexes)
        coord_names = self._coord_names & variables.keys() | coord_vars.keys()
        return self._replace_with_new_dims(variables, coord_names, indexes=indexes)

    def sel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        method: str | None = None,
        tolerance: int | float | Iterable[int | float] | None = None,
        drop: bool = False,
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new dataset with each array indexed by tick labels
        along the specified dimension(s).

        In contrast to `Dataset.isel`, indexers for this method should use
        labels instead of integers.

        Under the hood, this method is powered by using pandas's powerful Index
        objects. This makes label based indexing essentially just as fast as
        using integer indexing.

        It also means this method uses pandas's (well documented) logic for
        indexing. This means you can use string shortcuts for datetime indexes
        (e.g., '2000-01' to select all values in January 2000). It also means
        that slices are treated as inclusive of both the start and stop values,
        unlike normal Python indexing.

        Parameters
        ----------
        indexers : dict, optional
            A dict with keys matching dimensions and values given
            by scalars, slices or arrays of tick labels. For dimensions with
            multi-index, the indexer may also be a dict-like object with keys
            matching index level names.
            If DataArrays are passed as indexers, xarray-style indexing will be
            carried out. See :ref:`indexing` for the details.
            One of indexers or indexers_kwargs must be provided.
        method : {None, "nearest", "pad", "ffill", "backfill", "bfill"}, optional
            Method to use for inexact matches:

            * None (default): only exact matches
            * pad / ffill: propagate last valid index value forward
            * backfill / bfill: propagate next valid index value backward
            * nearest: use nearest valid index value
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations must
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
        drop : bool, optional
            If ``drop=True``, drop coordinates variables in `indexers` instead
            of making them scalar.
        **indexers_kwargs : {dim: indexer, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Returns
        -------
        obj : Dataset
            A new Dataset with the same contents as this dataset, except each
            variable and dimension is indexed by the appropriate indexers.
            If indexer DataArrays have coordinates that do not conflict with
            this object, then these coordinates will be attached.
            In general, each array's data will be a view of the array's data
            in this dataset, unless vectorized indexing was triggered by using
            an array indexer, in which case the data will be a copy.

        See Also
        --------
        :func:`Dataset.isel <Dataset.isel>`
        :func:`DataArray.sel <DataArray.sel>`

        :doc:`xarray-tutorial:intermediate/indexing/indexing`
            Tutorial material on indexing with Xarray objects

        :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
            Tutorial material on basics of indexing

        """
        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "sel")
        query_results = map_index_queries(
            self, indexers=indexers, method=method, tolerance=tolerance
        )

        if drop:
            no_scalar_variables = {}
            for k, v in query_results.variables.items():
                if v.dims:
                    no_scalar_variables[k] = v
                elif k in self._coord_names:
                    query_results.drop_coords.append(k)
            query_results.variables = no_scalar_variables

        result = self.isel(indexers=query_results.dim_indexers, drop=drop)
        return result._overwrite_indexes(*query_results.as_tuple()[1:])

    def _shuffle(self, dim, *, indices: GroupIndices, chunks: T_Chunks) -> Self:
        # Shuffling is only different from `isel` for chunked arrays.
        # Extract them out, and treat them specially. The rest, we route through isel.
        # This makes it easy to ensure correct handling of indexes.
        is_chunked = {
            name: var
            for name, var in self._variables.items()
            if is_chunked_array(var._data)
        }
        subset = self[[name for name in self._variables if name not in is_chunked]]

        no_slices: list[list[int]] = [
            (
                list(range(*idx.indices(self.sizes[dim])))
                if isinstance(idx, slice)
                else idx
            )
            for idx in indices
        ]
        no_slices = [idx for idx in no_slices if idx]

        shuffled = (
            subset
            if dim not in subset.dims
            else subset.isel({dim: np.concatenate(no_slices)})
        )
        for name, var in is_chunked.items():
            shuffled[name] = var._shuffle(
                indices=no_slices,
                dim=dim,
                chunks=chunks,
            )
        return shuffled

    def head(
        self,
        indexers: Mapping[Any, int] | int | None = None,
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new dataset with the first `n` values of each array
        for the specified dimension(s).

        Parameters
        ----------
        indexers : dict or int, default: 5
            A dict with keys matching dimensions and integer values `n`
            or a single integer `n` applied over all dimensions.
            One of indexers or indexers_kwargs must be provided.
        **indexers_kwargs : {dim: n, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Examples
        --------
        >>> dates = pd.date_range(start="2023-01-01", periods=5)
        >>> pageviews = [1200, 1500, 900, 1800, 2000]
        >>> visitors = [800, 1000, 600, 1200, 1500]
        >>> dataset = xr.Dataset(
        ...     {
        ...         "pageviews": (("date"), pageviews),
        ...         "visitors": (("date"), visitors),
        ...     },
        ...     coords={"date": dates},
        ... )
        >>> busiest_days = dataset.sortby("pageviews", ascending=False)
        >>> busiest_days.head()
        <xarray.Dataset> Size: 120B
        Dimensions:    (date: 5)
        Coordinates:
          * date       (date) datetime64[ns] 40B 2023-01-05 2023-01-04 ... 2023-01-03
        Data variables:
            pageviews  (date) int64 40B 2000 1800 1500 1200 900
            visitors   (date) int64 40B 1500 1200 1000 800 600

        # Retrieve the 3 most busiest days in terms of pageviews

        >>> busiest_days.head(3)
        <xarray.Dataset> Size: 72B
        Dimensions:    (date: 3)
        Coordinates:
          * date       (date) datetime64[ns] 24B 2023-01-05 2023-01-04 2023-01-02
        Data variables:
            pageviews  (date) int64 24B 2000 1800 1500
            visitors   (date) int64 24B 1500 1200 1000

        # Using a dictionary to specify the number of elements for specific dimensions

        >>> busiest_days.head({"date": 3})
        <xarray.Dataset> Size: 72B
        Dimensions:    (date: 3)
        Coordinates:
          * date       (date) datetime64[ns] 24B 2023-01-05 2023-01-04 2023-01-02
        Data variables:
            pageviews  (date) int64 24B 2000 1800 1500
            visitors   (date) int64 24B 1500 1200 1000

        See Also
        --------
        Dataset.tail
        Dataset.thin
        DataArray.head
        """
        if not indexers_kwargs:
            if indexers is None:
                indexers = 5
            if not isinstance(indexers, int) and not is_dict_like(indexers):
                raise TypeError("indexers must be either dict-like or a single integer")
        if isinstance(indexers, int):
            indexers = dict.fromkeys(self.dims, indexers)
        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "head")
        for k, v in indexers.items():
            if not isinstance(v, int):
                raise TypeError(
                    "expected integer type indexer for "
                    f"dimension {k!r}, found {type(v)!r}"
                )
            elif v < 0:
                raise ValueError(
                    "expected positive integer as indexer "
                    f"for dimension {k!r}, found {v}"
                )
        indexers_slices = {k: slice(val) for k, val in indexers.items()}
        return self.isel(indexers_slices)

    def tail(
        self,
        indexers: Mapping[Any, int] | int | None = None,
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new dataset with the last `n` values of each array
        for the specified dimension(s).

        Parameters
        ----------
        indexers : dict or int, default: 5
            A dict with keys matching dimensions and integer values `n`
            or a single integer `n` applied over all dimensions.
            One of indexers or indexers_kwargs must be provided.
        **indexers_kwargs : {dim: n, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Examples
        --------
        >>> activity_names = ["Walking", "Running", "Cycling", "Swimming", "Yoga"]
        >>> durations = [30, 45, 60, 45, 60]  # in minutes
        >>> energies = [150, 300, 250, 400, 100]  # in calories
        >>> dataset = xr.Dataset(
        ...     {
        ...         "duration": (["activity"], durations),
        ...         "energy_expenditure": (["activity"], energies),
        ...     },
        ...     coords={"activity": activity_names},
        ... )
        >>> sorted_dataset = dataset.sortby("energy_expenditure", ascending=False)
        >>> sorted_dataset
        <xarray.Dataset> Size: 240B
        Dimensions:             (activity: 5)
        Coordinates:
          * activity            (activity) <U8 160B 'Swimming' 'Running' ... 'Yoga'
        Data variables:
            duration            (activity) int64 40B 45 45 60 30 60
            energy_expenditure  (activity) int64 40B 400 300 250 150 100

        # Activities with the least energy expenditures using tail()

        >>> sorted_dataset.tail(3)
        <xarray.Dataset> Size: 144B
        Dimensions:             (activity: 3)
        Coordinates:
          * activity            (activity) <U8 96B 'Cycling' 'Walking' 'Yoga'
        Data variables:
            duration            (activity) int64 24B 60 30 60
            energy_expenditure  (activity) int64 24B 250 150 100

        >>> sorted_dataset.tail({"activity": 3})
        <xarray.Dataset> Size: 144B
        Dimensions:             (activity: 3)
        Coordinates:
          * activity            (activity) <U8 96B 'Cycling' 'Walking' 'Yoga'
        Data variables:
            duration            (activity) int64 24B 60 30 60
            energy_expenditure  (activity) int64 24B 250 150 100

        See Also
        --------
        Dataset.head
        Dataset.thin
        DataArray.tail
        """
        if not indexers_kwargs:
            if indexers is None:
                indexers = 5
            if not isinstance(indexers, int) and not is_dict_like(indexers):
                raise TypeError("indexers must be either dict-like or a single integer")
        if isinstance(indexers, int):
            indexers = dict.fromkeys(self.dims, indexers)
        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "tail")
        for k, v in indexers.items():
            if not isinstance(v, int):
                raise TypeError(
                    "expected integer type indexer for "
                    f"dimension {k!r}, found {type(v)!r}"
                )
            elif v < 0:
                raise ValueError(
                    "expected positive integer as indexer "
                    f"for dimension {k!r}, found {v}"
                )
        indexers_slices = {
            k: slice(-val, None) if val != 0 else slice(val)
            for k, val in indexers.items()
        }
        return self.isel(indexers_slices)

    def thin(
        self,
        indexers: Mapping[Any, int] | int | None = None,
        **indexers_kwargs: Any,
    ) -> Self:
        """Returns a new dataset with each array indexed along every `n`-th
        value for the specified dimension(s)

        Parameters
        ----------
        indexers : dict or int
            A dict with keys matching dimensions and integer values `n`
            or a single integer `n` applied over all dimensions.
            One of indexers or indexers_kwargs must be provided.
        **indexers_kwargs : {dim: n, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Examples
        --------
        >>> x_arr = np.arange(0, 26)
        >>> x_arr
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
               17, 18, 19, 20, 21, 22, 23, 24, 25])
        >>> x = xr.DataArray(
        ...     np.reshape(x_arr, (2, 13)),
        ...     dims=("x", "y"),
        ...     coords={"x": [0, 1], "y": np.arange(0, 13)},
        ... )
        >>> x_ds = xr.Dataset({"foo": x})
        >>> x_ds
        <xarray.Dataset> Size: 328B
        Dimensions:  (x: 2, y: 13)
        Coordinates:
          * x        (x) int64 16B 0 1
          * y        (y) int64 104B 0 1 2 3 4 5 6 7 8 9 10 11 12
        Data variables:
            foo      (x, y) int64 208B 0 1 2 3 4 5 6 7 8 ... 17 18 19 20 21 22 23 24 25

        >>> x_ds.thin(3)
        <xarray.Dataset> Size: 88B
        Dimensions:  (x: 1, y: 5)
        Coordinates:
          * x        (x) int64 8B 0
          * y        (y) int64 40B 0 3 6 9 12
        Data variables:
            foo      (x, y) int64 40B 0 3 6 9 12
        >>> x.thin({"x": 2, "y": 5})
        <xarray.DataArray (x: 1, y: 3)> Size: 24B
        array([[ 0,  5, 10]])
        Coordinates:
          * x        (x) int64 8B 0
          * y        (y) int64 24B 0 5 10

        See Also
        --------
        Dataset.head
        Dataset.tail
        DataArray.thin
        """
        if (
            not indexers_kwargs
            and not isinstance(indexers, int)
            and not is_dict_like(indexers)
        ):
            raise TypeError("indexers must be either dict-like or a single integer")
        if isinstance(indexers, int):
            indexers = dict.fromkeys(self.dims, indexers)
        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "thin")
        for k, v in indexers.items():
            if not isinstance(v, int):
                raise TypeError(
                    "expected integer type indexer for "
                    f"dimension {k!r}, found {type(v)!r}"
                )
            elif v < 0:
                raise ValueError(
                    "expected positive integer as indexer "
                    f"for dimension {k!r}, found {v}"
                )
            elif v == 0:
                raise ValueError("step cannot be zero")
        indexers_slices = {k: slice(None, None, val) for k, val in indexers.items()}
        return self.isel(indexers_slices)

    def broadcast_like(
        self,
        other: T_DataArrayOrSet,
        exclude: Iterable[Hashable] | None = None,
    ) -> Self:
        """Broadcast this DataArray against another Dataset or DataArray.
        This is equivalent to xr.broadcast(other, self)[1]

        Parameters
        ----------
        other : Dataset or DataArray
            Object against which to broadcast this array.
        exclude : iterable of hashable, optional
            Dimensions that must not be broadcasted

        """
        if exclude is None:
            exclude = set()
        else:
            exclude = set(exclude)
        args = align(other, self, join="outer", copy=False, exclude=exclude)

        dims_map, common_coords = _get_broadcast_dims_map_common_coords(args, exclude)

        return _broadcast_helper(args[1], exclude, dims_map, common_coords)

    def _reindex_callback(
        self,
        aligner: alignment.Aligner,
        dim_pos_indexers: dict[Hashable, Any],
        variables: dict[Hashable, Variable],
        indexes: dict[Hashable, Index],
        fill_value: Any,
        exclude_dims: frozenset[Hashable],
        exclude_vars: frozenset[Hashable],
    ) -> Self:
        """Callback called from ``Aligner`` to create a new reindexed Dataset."""

        new_variables = variables.copy()
        new_indexes = indexes.copy()

        # re-assign variable metadata
        for name, new_var in new_variables.items():
            var = self._variables.get(name)
            if var is not None:
                new_var.attrs = var.attrs
                new_var.encoding = var.encoding

        # pass through indexes from excluded dimensions
        # no extra check needed for multi-coordinate indexes, potential conflicts
        # should already have been detected when aligning the indexes
        for name, idx in self._indexes.items():
            var = self._variables[name]
            if set(var.dims) <= exclude_dims:
                new_indexes[name] = idx
                new_variables[name] = var

        if not dim_pos_indexers:
            # fast path for no reindexing necessary
            if set(new_indexes) - set(self._indexes):
                # this only adds new indexes and their coordinate variables
                reindexed = self._overwrite_indexes(new_indexes, new_variables)
            else:
                reindexed = self.copy(deep=aligner.copy)
        else:
            to_reindex = {
                k: v
                for k, v in self.variables.items()
                if k not in variables and k not in exclude_vars
            }
            reindexed_vars = alignment.reindex_variables(
                to_reindex,
                dim_pos_indexers,
                copy=aligner.copy,
                fill_value=fill_value,
                sparse=aligner.sparse,
            )
            new_variables.update(reindexed_vars)
            new_coord_names = self._coord_names | set(new_indexes)
            reindexed = self._replace_with_new_dims(
                new_variables, new_coord_names, indexes=new_indexes
            )

        reindexed.encoding = self.encoding

        return reindexed

    def reindex_like(
        self,
        other: T_Xarray,
        method: ReindexMethodOptions = None,
        tolerance: float | Iterable[float] | str | None = None,
        copy: bool = True,
        fill_value: Any = xrdtypes.NA,
    ) -> Self:
        """
        Conform this object onto the indexes of another object, for indexes which the
        objects share. Missing values are filled with ``fill_value``. The default fill
        value is NaN.

        Parameters
        ----------
        other : Dataset or DataArray
            Object with an 'indexes' attribute giving a mapping from dimension
            names to pandas.Index objects, which provides coordinates upon
            which to index the variables in this dataset. The indexes on this
            other object need not be the same as the indexes on this
            dataset. Any mismatched index values will be filled in with
            NaN, and any mismatched dimension names will simply be ignored.
        method : {None, "nearest", "pad", "ffill", "backfill", "bfill", None}, optional
            Method to use for filling index values from other not found in this
            dataset:

            - None (default): don't fill gaps
            - "pad" / "ffill": propagate last valid index value forward
            - "backfill" / "bfill": propagate next valid index value backward
            - "nearest": use nearest valid index value

        tolerance : float | Iterable[float] | str | None, default: None
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations must
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
            Tolerance may be a scalar value, which applies the same tolerance
            to all values, or list-like, which applies variable tolerance per
            element. List-like must be the same size as the index and its dtype
            must exactly match the index’s type.
        copy : bool, default: True
            If ``copy=True``, data in the return value is always copied. If
            ``copy=False`` and reindexing is unnecessary, or can be performed
            with only slice operations, then the output may share memory with
            the input. In either case, a new xarray object is always returned.
        fill_value : scalar or dict-like, optional
            Value to use for newly missing values. If a dict-like maps
            variable names to fill values.

        Returns
        -------
        reindexed : Dataset
            Another dataset, with this dataset's data but coordinates from the
            other object.

        See Also
        --------
        Dataset.reindex
        DataArray.reindex_like
        align

        """
        return alignment.reindex_like(
            self,
            other=other,
            method=method,
            tolerance=tolerance,
            copy=copy,
            fill_value=fill_value,
        )

    def reindex(
        self,
        indexers: Mapping[Any, Any] | None = None,
        method: ReindexMethodOptions = None,
        tolerance: float | Iterable[float] | str | None = None,
        copy: bool = True,
        fill_value: Any = xrdtypes.NA,
        **indexers_kwargs: Any,
    ) -> Self:
        """Conform this object onto a new set of indexes, filling in
        missing values with ``fill_value``. The default fill value is NaN.

        Parameters
        ----------
        indexers : dict, optional
            Dictionary with keys given by dimension names and values given by
            arrays of coordinates tick labels. Any mismatched coordinate
            values will be filled in with NaN, and any mismatched dimension
            names will simply be ignored.
            One of indexers or indexers_kwargs must be provided.
        method : {None, "nearest", "pad", "ffill", "backfill", "bfill", None}, optional
            Method to use for filling index values in ``indexers`` not found in
            this dataset:

            - None (default): don't fill gaps
            - "pad" / "ffill": propagate last valid index value forward
            - "backfill" / "bfill": propagate next valid index value backward
            - "nearest": use nearest valid index value

        tolerance : float | Iterable[float] | str | None, default: None
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations must
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
            Tolerance may be a scalar value, which applies the same tolerance
            to all values, or list-like, which applies variable tolerance per
            element. List-like must be the same size as the index and its dtype
            must exactly match the index’s type.
        copy : bool, default: True
            If ``copy=True``, data in the return value is always copied. If
            ``copy=False`` and reindexing is unnecessary, or can be performed
            with only slice operations, then the output may share memory with
            the input. In either case, a new xarray object is always returned.
        fill_value : scalar or dict-like, optional
            Value to use for newly missing values. If a dict-like,
            maps variable names (including coordinates) to fill values.
        sparse : bool, default: False
            use sparse-array.
        **indexers_kwargs : {dim: indexer, ...}, optional
            Keyword arguments in the same form as ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Returns
        -------
        reindexed : Dataset
            Another dataset, with this dataset's data but replaced coordinates.

        See Also
        --------
        Dataset.reindex_like
        align
        pandas.Index.get_indexer

        Examples
        --------
        Create a dataset with some fictional data.

        >>> x = xr.Dataset(
        ...     {
        ...         "temperature": ("station", 20 * np.random.rand(4)),
        ...         "pressure": ("station", 500 * np.random.rand(4)),
        ...     },
        ...     coords={"station": ["boston", "nyc", "seattle", "denver"]},
        ... )
        >>> x
        <xarray.Dataset> Size: 176B
        Dimensions:      (station: 4)
        Coordinates:
          * station      (station) <U7 112B 'boston' 'nyc' 'seattle' 'denver'
        Data variables:
            temperature  (station) float64 32B 10.98 14.3 12.06 10.9
            pressure     (station) float64 32B 211.8 322.9 218.8 445.9
        >>> x.indexes
        Indexes:
            station  Index(['boston', 'nyc', 'seattle', 'denver'], dtype='object', name='station')

        Create a new index and reindex the dataset. By default values in the new index that
        do not have corresponding records in the dataset are assigned `NaN`.

        >>> new_index = ["boston", "austin", "seattle", "lincoln"]
        >>> x.reindex({"station": new_index})
        <xarray.Dataset> Size: 176B
        Dimensions:      (station: 4)
        Coordinates:
          * station      (station) <U7 112B 'boston' 'austin' 'seattle' 'lincoln'
        Data variables:
            temperature  (station) float64 32B 10.98 nan 12.06 nan
            pressure     (station) float64 32B 211.8 nan 218.8 nan

        We can fill in the missing values by passing a value to the keyword `fill_value`.

        >>> x.reindex({"station": new_index}, fill_value=0)
        <xarray.Dataset> Size: 176B
        Dimensions:      (station: 4)
        Coordinates:
          * station      (station) <U7 112B 'boston' 'austin' 'seattle' 'lincoln'
        Data variables:
            temperature  (station) float64 32B 10.98 0.0 12.06 0.0
            pressure     (station) float64 32B 211.8 0.0 218.8 0.0

        We can also use different fill values for each variable.

        >>> x.reindex(
        ...     {"station": new_index}, fill_value={"temperature": 0, "pressure": 100}
        ... )
        <xarray.Dataset> Size: 176B
        Dimensions:      (station: 4)
        Coordinates:
          * station      (station) <U7 112B 'boston' 'austin' 'seattle' 'lincoln'
        Data variables:
            temperature  (station) float64 32B 10.98 0.0 12.06 0.0
            pressure     (station) float64 32B 211.8 100.0 218.8 100.0

        Because the index is not monotonically increasing or decreasing, we cannot use arguments
        to the keyword method to fill the `NaN` values.

        >>> x.reindex({"station": new_index}, method="nearest")
        Traceback (most recent call last):
        ...
            raise ValueError('index must be monotonic increasing or decreasing')
        ValueError: index must be monotonic increasing or decreasing

        To further illustrate the filling functionality in reindex, we will create a
        dataset with a monotonically increasing index (for example, a sequence of dates).

        >>> x2 = xr.Dataset(
        ...     {
        ...         "temperature": (
        ...             "time",
        ...             [15.57, 12.77, np.nan, 0.3081, 16.59, 15.12],
        ...         ),
        ...         "pressure": ("time", 500 * np.random.rand(6)),
        ...     },
        ...     coords={"time": pd.date_range("01/01/2019", periods=6, freq="D")},
        ... )
        >>> x2
        <xarray.Dataset> Size: 144B
        Dimensions:      (time: 6)
        Coordinates:
          * time         (time) datetime64[ns] 48B 2019-01-01 2019-01-02 ... 2019-01-06
        Data variables:
            temperature  (time) float64 48B 15.57 12.77 nan 0.3081 16.59 15.12
            pressure     (time) float64 48B 481.8 191.7 395.9 264.4 284.0 462.8

        Suppose we decide to expand the dataset to cover a wider date range.

        >>> time_index2 = pd.date_range("12/29/2018", periods=10, freq="D")
        >>> x2.reindex({"time": time_index2})
        <xarray.Dataset> Size: 240B
        Dimensions:      (time: 10)
        Coordinates:
          * time         (time) datetime64[ns] 80B 2018-12-29 2018-12-30 ... 2019-01-07
        Data variables:
            temperature  (time) float64 80B nan nan nan 15.57 ... 0.3081 16.59 15.12 nan
            pressure     (time) float64 80B nan nan nan 481.8 ... 264.4 284.0 462.8 nan

        The index entries that did not have a value in the original data frame (for example, `2018-12-29`)
        are by default filled with NaN. If desired, we can fill in the missing values using one of several options.

        For example, to back-propagate the last valid value to fill the `NaN` values,
        pass `bfill` as an argument to the `method` keyword.

        >>> x3 = x2.reindex({"time": time_index2}, method="bfill")
        >>> x3
        <xarray.Dataset> Size: 240B
        Dimensions:      (time: 10)
        Coordinates:
          * time         (time) datetime64[ns] 80B 2018-12-29 2018-12-30 ... 2019-01-07
        Data variables:
            temperature  (time) float64 80B 15.57 15.57 15.57 15.57 ... 16.59 15.12 nan
            pressure     (time) float64 80B 481.8 481.8 481.8 481.8 ... 284.0 462.8 nan

        Please note that the `NaN` value present in the original dataset (at index value `2019-01-03`)
        will not be filled by any of the value propagation schemes.

        >>> x2.where(x2.temperature.isnull(), drop=True)
        <xarray.Dataset> Size: 24B
        Dimensions:      (time: 1)
        Coordinates:
          * time         (time) datetime64[ns] 8B 2019-01-03
        Data variables:
            temperature  (time) float64 8B nan
            pressure     (time) float64 8B 395.9
        >>> x3.where(x3.temperature.isnull(), drop=True)
        <xarray.Dataset> Size: 48B
        Dimensions:      (time: 2)
        Coordinates:
          * time         (time) datetime64[ns] 16B 2019-01-03 2019-01-07
        Data variables:
            temperature  (time) float64 16B nan nan
            pressure     (time) float64 16B 395.9 nan

        This is because filling while reindexing does not look at dataset values, but only compares
        the original and desired indexes. If you do want to fill in the `NaN` values present in the
        original dataset, use the :py:meth:`~Dataset.fillna()` method.

        """
        indexers = utils.either_dict_or_kwargs(indexers, indexers_kwargs, "reindex")
        return alignment.reindex(
            self,
            indexers=indexers,
            method=method,
            tolerance=tolerance,
            copy=copy,
            fill_value=fill_value,
        )

    def _reindex(
        self,
        indexers: Mapping[Any, Any] | None = None,
        method: str | None = None,
        tolerance: int | float | Iterable[int | float] | None = None,
        copy: bool = True,
        fill_value: Any = xrdtypes.NA,
        sparse: bool = False,
        **indexers_kwargs: Any,
    ) -> Self:
        """
        Same as reindex but supports sparse option.
        """
        indexers = utils.either_dict_or_kwargs(indexers, indexers_kwargs, "reindex")
        return alignment.reindex(
            self,
            indexers=indexers,
            method=method,
            tolerance=tolerance,
            copy=copy,
            fill_value=fill_value,
            sparse=sparse,
        )

    def interp(
        self,
        coords: Mapping[Any, Any] | None = None,
        method: InterpOptions = "linear",
        assume_sorted: bool = False,
        kwargs: Mapping[str, Any] | None = None,
        method_non_numeric: str = "nearest",
        **coords_kwargs: Any,
    ) -> Self:
        """
        Interpolate a Dataset onto new coordinates.

        Performs univariate or multivariate interpolation of a Dataset onto new coordinates,
        utilizing either NumPy or SciPy interpolation routines.

        Out-of-range values are filled with NaN, unless specified otherwise via `kwargs` to the numpy/scipy interpolant.

        Parameters
        ----------
        coords : dict, optional
            Mapping from dimension names to the new coordinates.
            New coordinate can be a scalar, array-like or DataArray.
            If DataArrays are passed as new coordinates, their dimensions are
            used for the broadcasting. Missing values are skipped.
        method : { "linear", "nearest", "zero", "slinear", "quadratic", "cubic", \
            "quintic", "polynomial", "pchip", "barycentric", "krogh", "akima", "makima" }
            Interpolation method to use (see descriptions above).
        assume_sorted : bool, default: False
            If False, values of coordinates that are interpolated over can be
            in any order and they are sorted first. If True, interpolated
            coordinates are assumed to be an array of monotonically increasing
            values.
        kwargs : dict, optional
            Additional keyword arguments passed to the interpolator. Valid
            options and their behavior depend which interpolant is used.
        method_non_numeric : {"nearest", "pad", "ffill", "backfill", "bfill"}, optional
            Method for non-numeric types. Passed on to :py:meth:`Dataset.reindex`.
            ``"nearest"`` is used by default.
        **coords_kwargs : {dim: coordinate, ...}, optional
            The keyword arguments form of ``coords``.
            One of coords or coords_kwargs must be provided.


        Returns
        -------
        interpolated : Dataset
            New dataset on the new coordinates.

        Notes
        -----
        - SciPy is required for certain interpolation methods.
        - When interpolating along multiple dimensions with methods `linear` and `nearest`,
            the process attempts to decompose the interpolation into independent interpolations
            along one dimension at a time.
        - The specific interpolation method and dimensionality determine which
            interpolant is used:

            1. **Interpolation along one dimension of 1D data (`method='linear'`)**
                - Uses :py:func:`numpy.interp`, unless `fill_value='extrapolate'` is provided via `kwargs`.

            2. **Interpolation along one dimension of N-dimensional data (N ≥ 1)**
                - Methods {"linear", "nearest", "zero", "slinear", "quadratic", "cubic", "quintic", "polynomial"}
                    use :py:func:`scipy.interpolate.interp1d`, unless conditions permit the use of :py:func:`numpy.interp`
                    (as in the case of `method='linear'` for 1D data).
                - If `method='polynomial'`, the `order` keyword argument must also be provided.

            3. **Special interpolants for interpolation along one dimension of N-dimensional data (N ≥ 1)**
                - Depending on the `method`, the following interpolants from :py:class:`scipy.interpolate` are used:
                    - `"pchip"`: :py:class:`scipy.interpolate.PchipInterpolator`
                    - `"barycentric"`: :py:class:`scipy.interpolate.BarycentricInterpolator`
                    - `"krogh"`: :py:class:`scipy.interpolate.KroghInterpolator`
                    - `"akima"` or `"makima"`: :py:class:`scipy.interpolate.Akima1dInterpolator`
                        (`makima` is handled by passing the `makima` flag).

            4. **Interpolation along multiple dimensions of multi-dimensional data**
                - Uses :py:func:`scipy.interpolate.interpn` for methods {"linear", "nearest", "slinear",
                    "cubic", "quintic", "pchip"}.

        See Also
        --------
        :mod:`scipy.interpolate`

        :doc:`xarray-tutorial:fundamentals/02.2_manipulating_dimensions`
            Tutorial material on manipulating data resolution using :py:func:`~xarray.Dataset.interp`

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     data_vars={
        ...         "a": ("x", [5, 7, 4]),
        ...         "b": (
        ...             ("x", "y"),
        ...             [[1, 4, 2, 9], [2, 7, 6, np.nan], [6, np.nan, 5, 8]],
        ...         ),
        ...     },
        ...     coords={"x": [0, 1, 2], "y": [10, 12, 14, 16]},
        ... )
        >>> ds
        <xarray.Dataset> Size: 176B
        Dimensions:  (x: 3, y: 4)
        Coordinates:
          * x        (x) int64 24B 0 1 2
          * y        (y) int64 32B 10 12 14 16
        Data variables:
            a        (x) int64 24B 5 7 4
            b        (x, y) float64 96B 1.0 4.0 2.0 9.0 2.0 7.0 6.0 nan 6.0 nan 5.0 8.0

        1D interpolation with the default method (linear):

        >>> ds.interp(x=[0, 0.75, 1.25, 1.75])
        <xarray.Dataset> Size: 224B
        Dimensions:  (x: 4, y: 4)
        Coordinates:
          * y        (y) int64 32B 10 12 14 16
          * x        (x) float64 32B 0.0 0.75 1.25 1.75
        Data variables:
            a        (x) float64 32B 5.0 6.5 6.25 4.75
            b        (x, y) float64 128B 1.0 4.0 2.0 nan 1.75 ... nan 5.0 nan 5.25 nan

        1D interpolation with a different method:

        >>> ds.interp(x=[0, 0.75, 1.25, 1.75], method="nearest")
        <xarray.Dataset> Size: 224B
        Dimensions:  (x: 4, y: 4)
        Coordinates:
          * y        (y) int64 32B 10 12 14 16
          * x        (x) float64 32B 0.0 0.75 1.25 1.75
        Data variables:
            a        (x) float64 32B 5.0 7.0 7.0 4.0
            b        (x, y) float64 128B 1.0 4.0 2.0 9.0 2.0 7.0 ... nan 6.0 nan 5.0 8.0

        1D extrapolation:

        >>> ds.interp(
        ...     x=[1, 1.5, 2.5, 3.5],
        ...     method="linear",
        ...     kwargs={"fill_value": "extrapolate"},
        ... )
        <xarray.Dataset> Size: 224B
        Dimensions:  (x: 4, y: 4)
        Coordinates:
          * y        (y) int64 32B 10 12 14 16
          * x        (x) float64 32B 1.0 1.5 2.5 3.5
        Data variables:
            a        (x) float64 32B 7.0 5.5 2.5 -0.5
            b        (x, y) float64 128B 2.0 7.0 6.0 nan 4.0 ... nan 12.0 nan 3.5 nan

        2D interpolation:

        >>> ds.interp(x=[0, 0.75, 1.25, 1.75], y=[11, 13, 15], method="linear")
        <xarray.Dataset> Size: 184B
        Dimensions:  (x: 4, y: 3)
        Coordinates:
          * x        (x) float64 32B 0.0 0.75 1.25 1.75
          * y        (y) int64 24B 11 13 15
        Data variables:
            a        (x) float64 32B 5.0 6.5 6.25 4.75
            b        (x, y) float64 96B 2.5 3.0 nan 4.0 5.625 ... nan nan nan nan nan
        """
        from xarray.core import missing

        if kwargs is None:
            kwargs = {}

        coords = either_dict_or_kwargs(coords, coords_kwargs, "interp")
        indexers = dict(self._validate_interp_indexers(coords))
        obj = self if assume_sorted else self.sortby(list(coords))

        def maybe_variable(obj, k):
            # workaround to get variable for dimension without coordinate.
            try:
                return obj._variables[k]
            except KeyError:
                return as_variable((k, range(obj.sizes[k])))

        def _validate_interp_indexer(x, new_x):
            # In the case of datetimes, the restrictions placed on indexers
            # used with interp are stronger than those which are placed on
            # isel, so we need an additional check after _validate_indexers.
            if _contains_datetime_like_objects(
                x
            ) and not _contains_datetime_like_objects(new_x):
                raise TypeError(
                    "When interpolating over a datetime-like "
                    "coordinate, the coordinates to "
                    "interpolate to must be either datetime "
                    "strings or datetimes. "
                    f"Instead got\n{new_x}"
                )
            return x, new_x

        validated_indexers = {
            k: _validate_interp_indexer(maybe_variable(obj, k), v)
            for k, v in indexers.items()
        }

        # optimization: subset to coordinate range of the target index
        if method in ["linear", "nearest"]:
            for k, v in validated_indexers.items():
                obj, newidx = missing._localize(obj, {k: v})
                validated_indexers[k] = newidx[k]

        has_chunked_array = bool(
            any(is_chunked_array(v._data) for v in obj._variables.values())
        )
        if has_chunked_array:
            # optimization: create dask coordinate arrays once per Dataset
            # rather than once per Variable when dask.array.unify_chunks is called later
            # GH4739
            dask_indexers = {
                k: (index.to_base_variable().chunk(), dest.to_base_variable().chunk())
                for k, (index, dest) in validated_indexers.items()
            }

        variables: dict[Hashable, Variable] = {}
        reindex_vars: list[Hashable] = []
        for name, var in obj._variables.items():
            if name in indexers:
                continue

            use_indexers = (
                dask_indexers if is_duck_dask_array(var._data) else validated_indexers
            )

            dtype_kind = var.dtype.kind
            if dtype_kind in "uifc":
                # For normal number types do the interpolation:
                var_indexers = {k: v for k, v in use_indexers.items() if k in var.dims}
                variables[name] = missing.interp(var, var_indexers, method, **kwargs)
            elif dtype_kind in "ObU" and (use_indexers.keys() & var.dims):
                if all(var.sizes[d] == 1 for d in (use_indexers.keys() & var.dims)):
                    # Broadcastable, can be handled quickly without reindex:
                    to_broadcast = (var.squeeze(),) + tuple(
                        dest for _, dest in use_indexers.values()
                    )
                    variables[name] = broadcast_variables(*to_broadcast)[0].copy(
                        deep=True
                    )
                else:
                    # For types that we do not understand do stepwise
                    # interpolation to avoid modifying the elements.
                    # reindex the variable instead because it supports
                    # booleans and objects and retains the dtype but inside
                    # this loop there might be some duplicate code that slows it
                    # down, therefore collect these signals and run it later:
                    reindex_vars.append(name)
            elif all(d not in indexers for d in var.dims):
                # For anything else we can only keep variables if they
                # are not dependent on any coords that are being
                # interpolated along:
                variables[name] = var

        if reindex_vars and (
            reindex_indexers := {
                k: v for k, (_, v) in validated_indexers.items() if v.dims == (k,)
            }
        ):
            reindexed = alignment.reindex(
                obj[reindex_vars],
                indexers=reindex_indexers,
                method=method_non_numeric,
                exclude_vars=variables.keys(),
            )
            indexes = dict(reindexed._indexes)
            variables.update(reindexed.variables)
        else:
            # Get the indexes that are not being interpolated along
            indexes = {k: v for k, v in obj._indexes.items() if k not in indexers}

        # Get the coords that also exist in the variables:
        coord_names = obj._coord_names & variables.keys()
        selected = self._replace_with_new_dims(
            variables.copy(), coord_names, indexes=indexes
        )

        # Attach indexer as coordinate
        for k, v in indexers.items():
            assert isinstance(v, Variable)
            if v.dims == (k,):
                index = PandasIndex(v, k, coord_dtype=v.dtype)
                index_vars = index.create_variables({k: v})
                indexes[k] = index
                variables.update(index_vars)
            else:
                variables[k] = v

        # Extract coordinates from indexers
        coord_vars, new_indexes = selected._get_indexers_coords_and_indexes(coords)
        variables.update(coord_vars)
        indexes.update(new_indexes)

        coord_names = obj._coord_names & variables.keys() | coord_vars.keys()
        return self._replace_with_new_dims(variables, coord_names, indexes=indexes)

    def interp_like(
        self,
        other: T_Xarray,
        method: InterpOptions = "linear",
        assume_sorted: bool = False,
        kwargs: Mapping[str, Any] | None = None,
        method_non_numeric: str = "nearest",
    ) -> Self:
        """Interpolate this object onto the coordinates of another object.

        Performs univariate or multivariate interpolation of a Dataset onto new coordinates,
        utilizing either NumPy or SciPy interpolation routines.

        Out-of-range values are filled with NaN, unless specified otherwise via `kwargs` to the numpy/scipy interpolant.

        Parameters
        ----------
        other : Dataset or DataArray
            Object with an 'indexes' attribute giving a mapping from dimension
            names to an 1d array-like, which provides coordinates upon
            which to index the variables in this dataset. Missing values are skipped.
        method : { "linear", "nearest", "zero", "slinear", "quadratic", "cubic", \
            "quintic", "polynomial", "pchip", "barycentric", "krogh", "akima", "makima" }
            Interpolation method to use (see descriptions above).
        assume_sorted : bool, default: False
            If False, values of coordinates that are interpolated over can be
            in any order and they are sorted first. If True, interpolated
            coordinates are assumed to be an array of monotonically increasing
            values.
        kwargs : dict, optional
            Additional keyword arguments passed to the interpolator. Valid
            options and their behavior depend which interpolant is use
        method_non_numeric : {"nearest", "pad", "ffill", "backfill", "bfill"}, optional
            Method for non-numeric types. Passed on to :py:meth:`Dataset.reindex`.
            ``"nearest"`` is used by default.

        Returns
        -------
        interpolated : Dataset
            Another dataset by interpolating this dataset's data along the
            coordinates of the other object.

        Notes
        -----
        - scipy is required.
        - If the dataset has object-type coordinates, reindex is used for these
            coordinates instead of the interpolation.
        - When interpolating along multiple dimensions with methods `linear` and `nearest`,
            the process attempts to decompose the interpolation into independent interpolations
            along one dimension at a time.
        - The specific interpolation method and dimensionality determine which
            interpolant is used:

            1. **Interpolation along one dimension of 1D data (`method='linear'`)**
                - Uses :py:func:`numpy.interp`, unless `fill_value='extrapolate'` is provided via `kwargs`.

            2. **Interpolation along one dimension of N-dimensional data (N ≥ 1)**
                - Methods {"linear", "nearest", "zero", "slinear", "quadratic", "cubic", "quintic", "polynomial"}
                    use :py:func:`scipy.interpolate.interp1d`, unless conditions permit the use of :py:func:`numpy.interp`
                    (as in the case of `method='linear'` for 1D data).
                - If `method='polynomial'`, the `order` keyword argument must also be provided.

            3. **Special interpolants for interpolation along one dimension of N-dimensional data (N ≥ 1)**
                - Depending on the `method`, the following interpolants from :py:class:`scipy.interpolate` are used:
                    - `"pchip"`: :py:class:`scipy.interpolate.PchipInterpolator`
                    - `"barycentric"`: :py:class:`scipy.interpolate.BarycentricInterpolator`
                    - `"krogh"`: :py:class:`scipy.interpolate.KroghInterpolator`
                    - `"akima"` or `"makima"`: :py:class:`scipy.interpolate.Akima1dInterpolator`
                        (`makima` is handled by passing the `makima` flag).

            4. **Interpolation along multiple dimensions of multi-dimensional data**
                - Uses :py:func:`scipy.interpolate.interpn` for methods {"linear", "nearest", "slinear",
                    "cubic", "quintic", "pchip"}.

        See Also
        --------
        :func:`Dataset.interp`
        :func:`Dataset.reindex_like`
        :mod:`scipy.interpolate`
        """
        if kwargs is None:
            kwargs = {}

        # pick only dimension coordinates with a single index
        coords: dict[Hashable, Variable] = {}
        other_indexes = other.xindexes
        for dim in self.dims:
            other_dim_coords = other_indexes.get_all_coords(dim, errors="ignore")
            if len(other_dim_coords) == 1:
                coords[dim] = other_dim_coords[dim]

        numeric_coords: dict[Hashable, Variable] = {}
        object_coords: dict[Hashable, Variable] = {}
        for k, v in coords.items():
            if v.dtype.kind in "uifcMm":
                numeric_coords[k] = v
            else:
                object_coords[k] = v

        ds = self
        if object_coords:
            # We do not support interpolation along object coordinate.
            # reindex instead.
            ds = self.reindex(object_coords)
        return ds.interp(
            coords=numeric_coords,
            method=method,
            assume_sorted=assume_sorted,
            kwargs=kwargs,
            method_non_numeric=method_non_numeric,
        )

    # Helper methods for rename()
    def _rename_vars(
        self, name_dict, dims_dict
    ) -> tuple[dict[Hashable, Variable], set[Hashable]]:
        variables = {}
        coord_names = set()
        for k, v in self.variables.items():
            var = v.copy(deep=False)
            var.dims = tuple(dims_dict.get(dim, dim) for dim in v.dims)
            name = name_dict.get(k, k)
            if name in variables:
                raise ValueError(f"the new name {name!r} conflicts")
            variables[name] = var
            if k in self._coord_names:
                coord_names.add(name)
        return variables, coord_names

    def _rename_dims(self, name_dict: Mapping[Any, Hashable]) -> dict[Hashable, int]:
        return {name_dict.get(k, k): v for k, v in self.sizes.items()}

    def _rename_indexes(
        self, name_dict: Mapping[Any, Hashable], dims_dict: Mapping[Any, Hashable]
    ) -> tuple[dict[Hashable, Index], dict[Hashable, Variable]]:
        if not self._indexes:
            return {}, {}

        indexes = {}
        variables = {}

        for index, coord_names in self.xindexes.group_by_index():
            new_index = index.rename(name_dict, dims_dict)
            new_coord_names = [name_dict.get(k, k) for k in coord_names]
            indexes.update(dict.fromkeys(new_coord_names, new_index))
            new_index_vars = new_index.create_variables(
                {
                    new: self._variables[old]
                    for old, new in zip(coord_names, new_coord_names, strict=True)
                }
            )
            variables.update(new_index_vars)

        return indexes, variables

    def _rename_all(
        self, name_dict: Mapping[Any, Hashable], dims_dict: Mapping[Any, Hashable]
    ) -> tuple[
        dict[Hashable, Variable],
        set[Hashable],
        dict[Hashable, int],
        dict[Hashable, Index],
    ]:
        variables, coord_names = self._rename_vars(name_dict, dims_dict)
        dims = self._rename_dims(dims_dict)

        indexes, index_vars = self._rename_indexes(name_dict, dims_dict)
        variables = {k: index_vars.get(k, v) for k, v in variables.items()}

        return variables, coord_names, dims, indexes

    def _rename(
        self,
        name_dict: Mapping[Any, Hashable] | None = None,
        **names: Hashable,
    ) -> Self:
        """Also used internally by DataArray so that the warning (if any)
        is raised at the right stack level.
        """
        name_dict = either_dict_or_kwargs(name_dict, names, "rename")
        for k in name_dict.keys():
            if k not in self and k not in self.dims:
                raise ValueError(
                    f"cannot rename {k!r} because it is not a "
                    "variable or dimension in this dataset"
                )

            create_dim_coord = False
            new_k = name_dict[k]

            if k == new_k:
                continue  # Same name, nothing to do

            if k in self.dims and new_k in self._coord_names:
                coord_dims = self._variables[name_dict[k]].dims
                if coord_dims == (k,):
                    create_dim_coord = True
            elif k in self._coord_names and new_k in self.dims:
                coord_dims = self._variables[k].dims
                if coord_dims == (new_k,):
                    create_dim_coord = True

            if create_dim_coord:
                warnings.warn(
                    f"rename {k!r} to {name_dict[k]!r} does not create an index "
                    "anymore. Try using swap_dims instead or use set_index "
                    "after rename to create an indexed coordinate.",
                    UserWarning,
                    stacklevel=3,
                )

        variables, coord_names, dims, indexes = self._rename_all(
            name_dict=name_dict, dims_dict=name_dict
        )
        return self._replace(variables, coord_names, dims=dims, indexes=indexes)

    def rename(
        self,
        name_dict: Mapping[Any, Hashable] | None = None,
        **names: Hashable,
    ) -> Self:
        """Returns a new object with renamed variables, coordinates and dimensions.

        Parameters
        ----------
        name_dict : dict-like, optional
            Dictionary whose keys are current variable, coordinate or dimension names and
            whose values are the desired names.
        **names : optional
            Keyword form of ``name_dict``.
            One of name_dict or names must be provided.

        Returns
        -------
        renamed : Dataset
            Dataset with renamed variables, coordinates and dimensions.

        See Also
        --------
        Dataset.swap_dims
        Dataset.rename_vars
        Dataset.rename_dims
        DataArray.rename
        """
        return self._rename(name_dict=name_dict, **names)

    def rename_dims(
        self,
        dims_dict: Mapping[Any, Hashable] | None = None,
        **dims: Hashable,
    ) -> Self:
        """Returns a new object with renamed dimensions only.

        Parameters
        ----------
        dims_dict : dict-like, optional
            Dictionary whose keys are current dimension names and
            whose values are the desired names. The desired names must
            not be the name of an existing dimension or Variable in the Dataset.
        **dims : optional
            Keyword form of ``dims_dict``.
            One of dims_dict or dims must be provided.

        Returns
        -------
        renamed : Dataset
            Dataset with renamed dimensions.

        See Also
        --------
        Dataset.swap_dims
        Dataset.rename
        Dataset.rename_vars
        DataArray.rename
        """
        dims_dict = either_dict_or_kwargs(dims_dict, dims, "rename_dims")
        for k, v in dims_dict.items():
            if k not in self.dims:
                raise ValueError(
                    f"cannot rename {k!r} because it is not found "
                    f"in the dimensions of this dataset {tuple(self.dims)}"
                )
            if v in self.dims or v in self:
                raise ValueError(
                    f"Cannot rename {k} to {v} because {v} already exists. "
                    "Try using swap_dims instead."
                )

        variables, coord_names, sizes, indexes = self._rename_all(
            name_dict={}, dims_dict=dims_dict
        )
        return self._replace(variables, coord_names, dims=sizes, indexes=indexes)

    def rename_vars(
        self,
        name_dict: Mapping[Any, Hashable] | None = None,
        **names: Hashable,
    ) -> Self:
        """Returns a new object with renamed variables including coordinates

        Parameters
        ----------
        name_dict : dict-like, optional
            Dictionary whose keys are current variable or coordinate names and
            whose values are the desired names.
        **names : optional
            Keyword form of ``name_dict``.
            One of name_dict or names must be provided.

        Returns
        -------
        renamed : Dataset
            Dataset with renamed variables including coordinates

        See Also
        --------
        Dataset.swap_dims
        Dataset.rename
        Dataset.rename_dims
        DataArray.rename
        """
        name_dict = either_dict_or_kwargs(name_dict, names, "rename_vars")
        for k in name_dict:
            if k not in self:
                raise ValueError(
                    f"cannot rename {k!r} because it is not a "
                    "variable or coordinate in this dataset"
                )
        variables, coord_names, dims, indexes = self._rename_all(
            name_dict=name_dict, dims_dict={}
        )
        return self._replace(variables, coord_names, dims=dims, indexes=indexes)

    def swap_dims(
        self, dims_dict: Mapping[Any, Hashable] | None = None, **dims_kwargs
    ) -> Self:
        """Returns a new object with swapped dimensions.

        Parameters
        ----------
        dims_dict : dict-like
            Dictionary whose keys are current dimension names and whose values
            are new names.
        **dims_kwargs : {existing_dim: new_dim, ...}, optional
            The keyword arguments form of ``dims_dict``.
            One of dims_dict or dims_kwargs must be provided.

        Returns
        -------
        swapped : Dataset
            Dataset with swapped dimensions.

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     data_vars={"a": ("x", [5, 7]), "b": ("x", [0.1, 2.4])},
        ...     coords={"x": ["a", "b"], "y": ("x", [0, 1])},
        ... )
        >>> ds
        <xarray.Dataset> Size: 56B
        Dimensions:  (x: 2)
        Coordinates:
          * x        (x) <U1 8B 'a' 'b'
            y        (x) int64 16B 0 1
        Data variables:
            a        (x) int64 16B 5 7
            b        (x) float64 16B 0.1 2.4

        >>> ds.swap_dims({"x": "y"})
        <xarray.Dataset> Size: 56B
        Dimensions:  (y: 2)
        Coordinates:
            x        (y) <U1 8B 'a' 'b'
          * y        (y) int64 16B 0 1
        Data variables:
            a        (y) int64 16B 5 7
            b        (y) float64 16B 0.1 2.4

        >>> ds.swap_dims({"x": "z"})
        <xarray.Dataset> Size: 56B
        Dimensions:  (z: 2)
        Coordinates:
            x        (z) <U1 8B 'a' 'b'
            y        (z) int64 16B 0 1
        Dimensions without coordinates: z
        Data variables:
            a        (z) int64 16B 5 7
            b        (z) float64 16B 0.1 2.4

        See Also
        --------
        Dataset.rename
        DataArray.swap_dims
        """
        # TODO: deprecate this method in favor of a (less confusing)
        # rename_dims() method that only renames dimensions.

        dims_dict = either_dict_or_kwargs(dims_dict, dims_kwargs, "swap_dims")
        for current_name, new_name in dims_dict.items():
            if current_name not in self.dims:
                raise ValueError(
                    f"cannot swap from dimension {current_name!r} because it is "
                    f"not one of the dimensions of this dataset {tuple(self.dims)}"
                )
            if new_name in self.variables and self.variables[new_name].dims != (
                current_name,
            ):
                raise ValueError(
                    f"replacement dimension {new_name!r} is not a 1D "
                    f"variable along the old dimension {current_name!r}"
                )

        result_dims = {dims_dict.get(dim, dim) for dim in self.dims}

        coord_names = self._coord_names.copy()
        coord_names.update({dim for dim in dims_dict.values() if dim in self.variables})

        variables: dict[Hashable, Variable] = {}
        indexes: dict[Hashable, Index] = {}
        for current_name, current_variable in self.variables.items():
            dims = tuple(dims_dict.get(dim, dim) for dim in current_variable.dims)
            var: Variable
            if current_name in result_dims:
                var = current_variable.to_index_variable()
                var.dims = dims
                if current_name in self._indexes:
                    indexes[current_name] = self._indexes[current_name]
                    variables[current_name] = var
                else:
                    index, index_vars = create_default_index_implicit(var)
                    indexes.update(dict.fromkeys(index_vars, index))
                    variables.update(index_vars)
                    coord_names.update(index_vars)
            else:
                var = current_variable.to_base_variable()
                var.dims = dims
                variables[current_name] = var

        return self._replace_with_new_dims(variables, coord_names, indexes=indexes)

    def expand_dims(
        self,
        dim: Hashable | Sequence[Hashable] | Mapping[Any, Any] | None = None,
        axis: int | Sequence[int] | None = None,
        create_index_for_new_dim: bool = True,
        **dim_kwargs: Any,
    ) -> Self:
        """Return a new object with an additional axis (or axes) inserted at
        the corresponding position in the array shape.  The new object is a
        view into the underlying array, not a copy.

        If dim is already a scalar coordinate, it will be promoted to a 1D
        coordinate consisting of a single value.

        The automatic creation of indexes to back new 1D coordinate variables
        controlled by the create_index_for_new_dim kwarg.

        Parameters
        ----------
        dim : hashable, sequence of hashable, mapping, or None
            Dimensions to include on the new variable. If provided as hashable
            or sequence of hashable, then dimensions are inserted with length
            1. If provided as a mapping, then the keys are the new dimensions
            and the values are either integers (giving the length of the new
            dimensions) or array-like (giving the coordinates of the new
            dimensions).
        axis : int, sequence of int, or None, default: None
            Axis position(s) where new axis is to be inserted (position(s) on
            the result array). If a sequence of integers is passed,
            multiple axes are inserted. In this case, dim arguments should be
            same length list. If axis=None is passed, all the axes will be
            inserted to the start of the result array.
        create_index_for_new_dim : bool, default: True
            Whether to create new ``PandasIndex`` objects when the object being expanded contains scalar variables with names in ``dim``.
        **dim_kwargs : int or sequence or ndarray
            The keywords are arbitrary dimensions being inserted and the values
            are either the lengths of the new dims (if int is given), or their
            coordinates. Note, this is an alternative to passing a dict to the
            dim kwarg and will only be used if dim is None.

        Returns
        -------
        expanded : Dataset
            This object, but with additional dimension(s).

        Examples
        --------
        >>> dataset = xr.Dataset({"temperature": ([], 25.0)})
        >>> dataset
        <xarray.Dataset> Size: 8B
        Dimensions:      ()
        Data variables:
            temperature  float64 8B 25.0

        # Expand the dataset with a new dimension called "time"

        >>> dataset.expand_dims(dim="time")
        <xarray.Dataset> Size: 8B
        Dimensions:      (time: 1)
        Dimensions without coordinates: time
        Data variables:
            temperature  (time) float64 8B 25.0

        # 1D data

        >>> temperature_1d = xr.DataArray([25.0, 26.5, 24.8], dims="x")
        >>> dataset_1d = xr.Dataset({"temperature": temperature_1d})
        >>> dataset_1d
        <xarray.Dataset> Size: 24B
        Dimensions:      (x: 3)
        Dimensions without coordinates: x
        Data variables:
            temperature  (x) float64 24B 25.0 26.5 24.8

        # Expand the dataset with a new dimension called "time" using axis argument

        >>> dataset_1d.expand_dims(dim="time", axis=0)
        <xarray.Dataset> Size: 24B
        Dimensions:      (time: 1, x: 3)
        Dimensions without coordinates: time, x
        Data variables:
            temperature  (time, x) float64 24B 25.0 26.5 24.8

        # 2D data

        >>> temperature_2d = xr.DataArray(np.random.rand(3, 4), dims=("y", "x"))
        >>> dataset_2d = xr.Dataset({"temperature": temperature_2d})
        >>> dataset_2d
        <xarray.Dataset> Size: 96B
        Dimensions:      (y: 3, x: 4)
        Dimensions without coordinates: y, x
        Data variables:
            temperature  (y, x) float64 96B 0.5488 0.7152 0.6028 ... 0.7917 0.5289

        # Expand the dataset with a new dimension called "time" using axis argument

        >>> dataset_2d.expand_dims(dim="time", axis=2)
        <xarray.Dataset> Size: 96B
        Dimensions:      (y: 3, x: 4, time: 1)
        Dimensions without coordinates: y, x, time
        Data variables:
            temperature  (y, x, time) float64 96B 0.5488 0.7152 0.6028 ... 0.7917 0.5289

        # Expand a scalar variable along a new dimension of the same name with and without creating a new index

        >>> ds = xr.Dataset(coords={"x": 0})
        >>> ds
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Coordinates:
            x        int64 8B 0
        Data variables:
            *empty*

        >>> ds.expand_dims("x")
        <xarray.Dataset> Size: 8B
        Dimensions:  (x: 1)
        Coordinates:
          * x        (x) int64 8B 0
        Data variables:
            *empty*

        >>> ds.expand_dims("x").indexes
        Indexes:
            x        Index([0], dtype='int64', name='x')

        >>> ds.expand_dims("x", create_index_for_new_dim=False).indexes
        Indexes:
            *empty*

        See Also
        --------
        DataArray.expand_dims
        """
        if dim is None:
            pass
        elif isinstance(dim, Mapping):
            # We're later going to modify dim in place; don't tamper with
            # the input
            dim = dict(dim)
        elif isinstance(dim, int):
            raise TypeError(
                "dim should be hashable or sequence of hashables or mapping"
            )
        elif isinstance(dim, str) or not isinstance(dim, Sequence):
            dim = {dim: 1}
        elif isinstance(dim, Sequence):
            if len(dim) != len(set(dim)):
                raise ValueError("dims should not contain duplicate values.")
            dim = dict.fromkeys(dim, 1)

        dim = either_dict_or_kwargs(dim, dim_kwargs, "expand_dims")
        assert isinstance(dim, MutableMapping)

        if axis is None:
            axis = list(range(len(dim)))
        elif not isinstance(axis, Sequence):
            axis = [axis]

        if len(dim) != len(axis):
            raise ValueError("lengths of dim and axis should be identical.")
        for d in dim:
            if d in self.dims:
                raise ValueError(f"Dimension {d} already exists.")
            if d in self._variables and not utils.is_scalar(self._variables[d]):
                raise ValueError(f"{d} already exists as coordinate or variable name.")

        variables: dict[Hashable, Variable] = {}
        indexes: dict[Hashable, Index] = dict(self._indexes)
        coord_names = self._coord_names.copy()
        # If dim is a dict, then ensure that the values are either integers
        # or iterables.
        for k, v in dim.items():
            if hasattr(v, "__iter__"):
                # If the value for the new dimension is an iterable, then
                # save the coordinates to the variables dict, and set the
                # value within the dim dict to the length of the iterable
                # for later use.

                if create_index_for_new_dim:
                    index = PandasIndex(v, k)
                    indexes[k] = index
                    name_and_new_1d_var = index.create_variables()
                else:
                    name_and_new_1d_var = {k: Variable(data=v, dims=k)}
                variables.update(name_and_new_1d_var)
                coord_names.add(k)
                dim[k] = variables[k].size
            elif isinstance(v, int):
                pass  # Do nothing if the dimensions value is just an int
            else:
                raise TypeError(
                    f"The value of new dimension {k} must be an iterable or an int"
                )

        for k, v in self._variables.items():
            if k not in dim:
                if k in coord_names:  # Do not change coordinates
                    variables[k] = v
                else:
                    result_ndim = len(v.dims) + len(axis)
                    for a in axis:
                        if a < -result_ndim or result_ndim - 1 < a:
                            raise IndexError(
                                f"Axis {a} of variable {k} is out of bounds of the "
                                f"expanded dimension size {result_ndim}"
                            )

                    axis_pos = [a if a >= 0 else result_ndim + a for a in axis]
                    if len(axis_pos) != len(set(axis_pos)):
                        raise ValueError("axis should not contain duplicate values")
                    # We need to sort them to make sure `axis` equals to the
                    # axis positions of the result array.
                    zip_axis_dim = sorted(zip(axis_pos, dim.items(), strict=True))

                    all_dims = list(zip(v.dims, v.shape, strict=True))
                    for d, c in zip_axis_dim:
                        all_dims.insert(d, c)
                    variables[k] = v.set_dims(dict(all_dims))
            elif k not in variables:
                if k in coord_names and create_index_for_new_dim:
                    # If dims includes a label of a non-dimension coordinate,
                    # it will be promoted to a 1D coordinate with a single value.
                    index, index_vars = create_default_index_implicit(v.set_dims(k))
                    indexes[k] = index
                    variables.update(index_vars)
                else:
                    if create_index_for_new_dim:
                        warnings.warn(
                            f"No index created for dimension {k} because variable {k} is not a coordinate. "
                            f"To create an index for {k}, please first call `.set_coords('{k}')` on this object.",
                            UserWarning,
                            stacklevel=2,
                        )

                    # create 1D variable without creating a new index
                    new_1d_var = v.set_dims(k)
                    variables.update({k: new_1d_var})

        return self._replace_with_new_dims(
            variables, coord_names=coord_names, indexes=indexes
        )

    def set_index(
        self,
        indexes: Mapping[Any, Hashable | Sequence[Hashable]] | None = None,
        append: bool = False,
        **indexes_kwargs: Hashable | Sequence[Hashable],
    ) -> Self:
        """Set Dataset (multi-)indexes using one or more existing coordinates
        or variables.

        This legacy method is limited to pandas (multi-)indexes and
        1-dimensional "dimension" coordinates. See
        :py:meth:`~Dataset.set_xindex` for setting a pandas or a custom
        Xarray-compatible index from one or more arbitrary coordinates.

        Parameters
        ----------
        indexes : {dim: index, ...}
            Mapping from names matching dimensions and values given
            by (lists of) the names of existing coordinates or variables to set
            as new (multi-)index.
        append : bool, default: False
            If True, append the supplied index(es) to the existing index(es).
            Otherwise replace the existing index(es) (default).
        **indexes_kwargs : optional
            The keyword arguments form of ``indexes``.
            One of indexes or indexes_kwargs must be provided.

        Returns
        -------
        obj : Dataset
            Another dataset, with this dataset's data but replaced coordinates.

        Examples
        --------
        >>> arr = xr.DataArray(
        ...     data=np.ones((2, 3)),
        ...     dims=["x", "y"],
        ...     coords={"x": range(2), "y": range(3), "a": ("x", [3, 4])},
        ... )
        >>> ds = xr.Dataset({"v": arr})
        >>> ds
        <xarray.Dataset> Size: 104B
        Dimensions:  (x: 2, y: 3)
        Coordinates:
          * x        (x) int64 16B 0 1
          * y        (y) int64 24B 0 1 2
            a        (x) int64 16B 3 4
        Data variables:
            v        (x, y) float64 48B 1.0 1.0 1.0 1.0 1.0 1.0
        >>> ds.set_index(x="a")
        <xarray.Dataset> Size: 88B
        Dimensions:  (x: 2, y: 3)
        Coordinates:
          * x        (x) int64 16B 3 4
          * y        (y) int64 24B 0 1 2
        Data variables:
            v        (x, y) float64 48B 1.0 1.0 1.0 1.0 1.0 1.0

        See Also
        --------
        Dataset.reset_index
        Dataset.set_xindex
        Dataset.swap_dims
        """
        dim_coords = either_dict_or_kwargs(indexes, indexes_kwargs, "set_index")

        new_indexes: dict[Hashable, Index] = {}
        new_variables: dict[Hashable, Variable] = {}
        drop_indexes: set[Hashable] = set()
        drop_variables: set[Hashable] = set()
        replace_dims: dict[Hashable, Hashable] = {}
        all_var_names: set[Hashable] = set()

        for dim, _var_names in dim_coords.items():
            if isinstance(_var_names, str) or not isinstance(_var_names, Sequence):
                var_names = [_var_names]
            else:
                var_names = list(_var_names)

            invalid_vars = set(var_names) - set(self._variables)
            if invalid_vars:
                raise ValueError(
                    ", ".join([str(v) for v in invalid_vars])
                    + " variable(s) do not exist"
                )

            all_var_names.update(var_names)
            drop_variables.update(var_names)

            # drop any pre-existing index involved and its corresponding coordinates
            index_coord_names = self.xindexes.get_all_coords(dim, errors="ignore")
            all_index_coord_names = set(index_coord_names)
            for k in var_names:
                all_index_coord_names.update(
                    self.xindexes.get_all_coords(k, errors="ignore")
                )

            drop_indexes.update(all_index_coord_names)
            drop_variables.update(all_index_coord_names)

            if len(var_names) == 1 and (not append or dim not in self._indexes):
                var_name = var_names[0]
                var = self._variables[var_name]
                # an error with a better message will be raised for scalar variables
                # when creating the PandasIndex
                if var.ndim > 0 and var.dims != (dim,):
                    raise ValueError(
                        f"dimension mismatch: try setting an index for dimension {dim!r} with "
                        f"variable {var_name!r} that has dimensions {var.dims}"
                    )
                idx = PandasIndex.from_variables({dim: var}, options={})
                idx_vars = idx.create_variables({var_name: var})

                # trick to preserve coordinate order in this case
                if dim in self._coord_names:
                    drop_variables.remove(dim)
            else:
                if append:
                    current_variables = {
                        k: self._variables[k] for k in index_coord_names
                    }
                else:
                    current_variables = {}
                idx, idx_vars = PandasMultiIndex.from_variables_maybe_expand(
                    dim,
                    current_variables,
                    {k: self._variables[k] for k in var_names},
                )
                for n in idx.index.names:
                    replace_dims[n] = dim

            new_indexes.update(dict.fromkeys(idx_vars, idx))
            new_variables.update(idx_vars)

        # re-add deindexed coordinates (convert to base variables)
        for k in drop_variables:
            if (
                k not in new_variables
                and k not in all_var_names
                and k in self._coord_names
            ):
                new_variables[k] = self._variables[k].to_base_variable()

        indexes_: dict[Any, Index] = {
            k: v for k, v in self._indexes.items() if k not in drop_indexes
        }
        indexes_.update(new_indexes)

        variables = {
            k: v for k, v in self._variables.items() if k not in drop_variables
        }
        variables.update(new_variables)

        # update dimensions if necessary, GH: 3512
        for k, v in variables.items():
            if any(d in replace_dims for d in v.dims):
                new_dims = [replace_dims.get(d, d) for d in v.dims]
                variables[k] = v._replace(dims=new_dims)

        coord_names = self._coord_names - drop_variables | set(new_variables)

        return self._replace_with_new_dims(
            variables, coord_names=coord_names, indexes=indexes_
        )

    def reset_index(
        self,
        dims_or_levels: Hashable | Sequence[Hashable],
        *,
        drop: bool = False,
    ) -> Self:
        """Reset the specified index(es) or multi-index level(s).

        This legacy method is specific to pandas (multi-)indexes and
        1-dimensional "dimension" coordinates. See the more generic
        :py:meth:`~Dataset.drop_indexes` and :py:meth:`~Dataset.set_xindex`
        method to respectively drop and set pandas or custom indexes for
        arbitrary coordinates.

        Parameters
        ----------
        dims_or_levels : Hashable or Sequence of Hashable
            Name(s) of the dimension(s) and/or multi-index level(s) that will
            be reset.
        drop : bool, default: False
            If True, remove the specified indexes and/or multi-index levels
            instead of extracting them as new coordinates (default: False).

        Returns
        -------
        obj : Dataset
            Another dataset, with this dataset's data but replaced coordinates.

        See Also
        --------
        Dataset.set_index
        Dataset.set_xindex
        Dataset.drop_indexes
        """
        if isinstance(dims_or_levels, str) or not isinstance(dims_or_levels, Sequence):
            dims_or_levels = [dims_or_levels]

        invalid_coords = set(dims_or_levels) - set(self._indexes)
        if invalid_coords:
            raise ValueError(
                f"{tuple(invalid_coords)} are not coordinates with an index"
            )

        drop_indexes: set[Hashable] = set()
        drop_variables: set[Hashable] = set()
        seen: set[Index] = set()
        new_indexes: dict[Hashable, Index] = {}
        new_variables: dict[Hashable, Variable] = {}

        def drop_or_convert(var_names):
            if drop:
                drop_variables.update(var_names)
            else:
                base_vars = {
                    k: self._variables[k].to_base_variable() for k in var_names
                }
                new_variables.update(base_vars)

        for name in dims_or_levels:
            index = self._indexes[name]

            if index in seen:
                continue
            seen.add(index)

            idx_var_names = set(self.xindexes.get_all_coords(name))
            drop_indexes.update(idx_var_names)

            if isinstance(index, PandasMultiIndex):
                # special case for pd.MultiIndex
                level_names = index.index.names
                keep_level_vars = {
                    k: self._variables[k]
                    for k in level_names
                    if k not in dims_or_levels
                }

                if index.dim not in dims_or_levels and keep_level_vars:
                    # do not drop the multi-index completely
                    # instead replace it by a new (multi-)index with dropped level(s)
                    idx = index.keep_levels(keep_level_vars)
                    idx_vars = idx.create_variables(keep_level_vars)
                    new_indexes.update(dict.fromkeys(idx_vars, idx))
                    new_variables.update(idx_vars)
                    if not isinstance(idx, PandasMultiIndex):
                        # multi-index reduced to single index
                        # backward compatibility: unique level coordinate renamed to dimension
                        drop_variables.update(keep_level_vars)
                    drop_or_convert(
                        [k for k in level_names if k not in keep_level_vars]
                    )
                else:
                    # always drop the multi-index dimension variable
                    drop_variables.add(index.dim)
                    drop_or_convert(level_names)
            else:
                drop_or_convert(idx_var_names)

        indexes = {k: v for k, v in self._indexes.items() if k not in drop_indexes}
        indexes.update(new_indexes)

        variables = {
            k: v for k, v in self._variables.items() if k not in drop_variables
        }
        variables.update(new_variables)

        coord_names = self._coord_names - drop_variables

        return self._replace_with_new_dims(
            variables, coord_names=coord_names, indexes=indexes
        )

    def set_xindex(
        self,
        coord_names: str | Sequence[Hashable],
        index_cls: type[Index] | None = None,
        **options,
    ) -> Self:
        """Set a new, Xarray-compatible index from one or more existing
        coordinate(s).

        Parameters
        ----------
        coord_names : str or list
            Name(s) of the coordinate(s) used to build the index.
            If several names are given, their order matters.
        index_cls : subclass of :class:`~xarray.indexes.Index`, optional
            The type of index to create. By default, try setting
            a ``PandasIndex`` if ``len(coord_names) == 1``,
            otherwise a ``PandasMultiIndex``.
        **options
            Options passed to the index constructor.

        Returns
        -------
        obj : Dataset
            Another dataset, with this dataset's data and with a new index.

        """
        # the Sequence check is required for mypy
        if is_scalar(coord_names) or not isinstance(coord_names, Sequence):
            coord_names = [coord_names]

        if index_cls is None:
            if len(coord_names) == 1:
                index_cls = PandasIndex
            else:
                index_cls = PandasMultiIndex
        elif not issubclass(index_cls, Index):
            raise TypeError(f"{index_cls} is not a subclass of xarray.Index")

        invalid_coords = set(coord_names) - self._coord_names

        if invalid_coords:
            msg = ["invalid coordinate(s)"]
            no_vars = invalid_coords - set(self._variables)
            data_vars = invalid_coords - no_vars
            if no_vars:
                msg.append(f"those variables don't exist: {no_vars}")
            if data_vars:
                msg.append(
                    f"those variables are data variables: {data_vars}, use `set_coords` first"
                )
            raise ValueError("\n".join(msg))

        # we could be more clever here (e.g., drop-in index replacement if index
        # coordinates do not conflict), but let's not allow this for now
        indexed_coords = set(coord_names) & set(self._indexes)

        if indexed_coords:
            raise ValueError(
                f"those coordinates already have an index: {indexed_coords}"
            )

        coord_vars = {name: self._variables[name] for name in coord_names}

        index = index_cls.from_variables(coord_vars, options=options)

        new_coord_vars = index.create_variables(coord_vars)

        # special case for setting a pandas multi-index from level coordinates
        # TODO: remove it once we depreciate pandas multi-index dimension (tuple
        # elements) coordinate
        if isinstance(index, PandasMultiIndex):
            coord_names = [index.dim] + list(coord_names)

        # Check for extra variables that don't match the coordinate names
        extra_vars = set(new_coord_vars) - set(coord_names)
        if extra_vars:
            extra_vars_str = ", ".join(f"'{name}'" for name in extra_vars)
            coord_names_str = ", ".join(f"'{name}'" for name in coord_names)
            raise ValueError(
                f"The index created extra variables {extra_vars_str} that are not "
                f"in the list of coordinates {coord_names_str}. "
                f"Use a factory method pattern instead:\n"
                f"  index = {index_cls.__name__}.from_variables(ds, {list(coord_names)!r})\n"
                f"  coords = xr.Coordinates.from_xindex(index)\n"
                f"  ds = ds.assign_coords(coords)"
            )

        variables: dict[Hashable, Variable]
        indexes: dict[Hashable, Index]

        if len(coord_names) == 1:
            variables = self._variables.copy()
            indexes = self._indexes.copy()

            name = list(coord_names).pop()
            if name in new_coord_vars:
                variables[name] = new_coord_vars[name]
            indexes[name] = index
        else:
            # reorder variables and indexes so that coordinates having the same
            # index are next to each other
            variables = {}
            for name, var in self._variables.items():
                if name not in coord_names:
                    variables[name] = var

            indexes = {}
            for name, idx in self._indexes.items():
                if name not in coord_names:
                    indexes[name] = idx

            for name in coord_names:
                try:
                    variables[name] = new_coord_vars[name]
                except KeyError:
                    variables[name] = self._variables[name]
                indexes[name] = index

        return self._replace(
            variables=variables,
            coord_names=self._coord_names | set(coord_names),
            indexes=indexes,
        )

    def reorder_levels(
        self,
        dim_order: Mapping[Any, Sequence[int | Hashable]] | None = None,
        **dim_order_kwargs: Sequence[int | Hashable],
    ) -> Self:
        """Rearrange index levels using input order.

        Parameters
        ----------
        dim_order : dict-like of Hashable to Sequence of int or Hashable, optional
            Mapping from names matching dimensions and values given
            by lists representing new level orders. Every given dimension
            must have a multi-index.
        **dim_order_kwargs : Sequence of int or Hashable, optional
            The keyword arguments form of ``dim_order``.
            One of dim_order or dim_order_kwargs must be provided.

        Returns
        -------
        obj : Dataset
            Another dataset, with this dataset's data but replaced
            coordinates.
        """
        dim_order = either_dict_or_kwargs(dim_order, dim_order_kwargs, "reorder_levels")
        variables = self._variables.copy()
        indexes = dict(self._indexes)
        new_indexes: dict[Hashable, Index] = {}
        new_variables: dict[Hashable, IndexVariable] = {}

        for dim, order in dim_order.items():
            index = self._indexes[dim]

            if not isinstance(index, PandasMultiIndex):
                raise ValueError(f"coordinate {dim} has no MultiIndex")

            level_vars = {k: self._variables[k] for k in order}
            idx = index.reorder_levels(level_vars)
            idx_vars = idx.create_variables(level_vars)
            new_indexes.update(dict.fromkeys(idx_vars, idx))
            new_variables.update(idx_vars)

        indexes = {k: v for k, v in self._indexes.items() if k not in new_indexes}
        indexes.update(new_indexes)

        variables = {k: v for k, v in self._variables.items() if k not in new_variables}
        variables.update(new_variables)

        return self._replace(variables, indexes=indexes)

    def _get_stack_index(
        self,
        dim,
        multi=False,
        create_index=False,
    ) -> tuple[Index | None, dict[Hashable, Variable]]:
        """Used by stack and unstack to get one pandas (multi-)index among
        the indexed coordinates along dimension `dim`.

        If exactly one index is found, return it with its corresponding
        coordinate variables(s), otherwise return None and an empty dict.

        If `create_index=True`, create a new index if none is found or raise
        an error if multiple indexes are found.

        """
        stack_index: Index | None = None
        stack_coords: dict[Hashable, Variable] = {}

        for name, index in self._indexes.items():
            var = self._variables[name]
            if (
                var.ndim == 1
                and var.dims[0] == dim
                and (
                    # stack: must be a single coordinate index
                    (not multi and not self.xindexes.is_multi(name))
                    # unstack: must be an index that implements .unstack
                    or (multi and type(index).unstack is not Index.unstack)
                )
            ):
                if stack_index is not None and index is not stack_index:
                    # more than one index found, stop
                    if create_index:
                        raise ValueError(
                            f"cannot stack dimension {dim!r} with `create_index=True` "
                            "and with more than one index found along that dimension"
                        )
                    return None, {}
                stack_index = index
                stack_coords[name] = var

        if create_index and stack_index is None:
            if dim in self._variables:
                var = self._variables[dim]
            else:
                _, _, var = _get_virtual_variable(self._variables, dim, self.sizes)
            # dummy index (only `stack_coords` will be used to construct the multi-index)
            stack_index = PandasIndex([0], dim)
            stack_coords = {dim: var}

        return stack_index, stack_coords

    def _stack_once(
        self,
        dims: Sequence[Hashable | EllipsisType],
        new_dim: Hashable,
        index_cls: type[Index],
        create_index: bool | None = True,
    ) -> Self:
        if dims == ...:
            raise ValueError("Please use [...] for dims, rather than just ...")
        if ... in dims:
            dims = list(infix_dims(dims, self.dims))

        new_variables: dict[Hashable, Variable] = {}
        stacked_var_names: list[Hashable] = []
        drop_indexes: list[Hashable] = []

        for name, var in self.variables.items():
            if any(d in var.dims for d in dims):
                add_dims = [d for d in dims if d not in var.dims]
                vdims = list(var.dims) + add_dims
                shape = [self.sizes[d] for d in vdims]
                exp_var = var.set_dims(vdims, shape)
                stacked_var = exp_var.stack(**{new_dim: dims})
                new_variables[name] = stacked_var
                stacked_var_names.append(name)
            else:
                new_variables[name] = var.copy(deep=False)

        # drop indexes of stacked coordinates (if any)
        for name in stacked_var_names:
            drop_indexes += list(self.xindexes.get_all_coords(name, errors="ignore"))

        new_indexes = {}
        new_coord_names = set(self._coord_names)
        if create_index or create_index is None:
            product_vars: dict[Any, Variable] = {}
            for dim in dims:
                idx, idx_vars = self._get_stack_index(dim, create_index=create_index)
                if idx is not None:
                    product_vars.update(idx_vars)

            if len(product_vars) == len(dims):
                idx = index_cls.stack(product_vars, new_dim)
                new_indexes[new_dim] = idx
                new_indexes.update(dict.fromkeys(product_vars, idx))
                idx_vars = idx.create_variables(product_vars)
                # keep consistent multi-index coordinate order
                for k in idx_vars:
                    new_variables.pop(k, None)
                new_variables.update(idx_vars)
                new_coord_names.update(idx_vars)

        indexes = {k: v for k, v in self._indexes.items() if k not in drop_indexes}
        indexes.update(new_indexes)

        return self._replace_with_new_dims(
            new_variables, coord_names=new_coord_names, indexes=indexes
        )

    @partial(deprecate_dims, old_name="dimensions")
    def stack(
        self,
        dim: Mapping[Any, Sequence[Hashable | EllipsisType]] | None = None,
        create_index: bool | None = True,
        index_cls: type[Index] = PandasMultiIndex,
        **dim_kwargs: Sequence[Hashable | EllipsisType],
    ) -> Self:
        """
        Stack any number of existing dimensions into a single new dimension.

        New dimensions will be added at the end, and by default the corresponding
        coordinate variables will be combined into a MultiIndex.

        Parameters
        ----------
        dim : mapping of hashable to sequence of hashable
            Mapping of the form `new_name=(dim1, dim2, ...)`. Names of new
            dimensions, and the existing dimensions that they replace. An
            ellipsis (`...`) will be replaced by all unlisted dimensions.
            Passing a list containing an ellipsis (`stacked_dim=[...]`) will stack over
            all dimensions.
        create_index : bool or None, default: True

            - True: create a multi-index for each of the stacked dimensions.
            - False: don't create any index.
            - None. create a multi-index only if exactly one single (1-d) coordinate
              index is found for every dimension to stack.

        index_cls: Index-class, default: PandasMultiIndex
            Can be used to pass a custom multi-index type (must be an Xarray index that
            implements `.stack()`). By default, a pandas multi-index wrapper is used.
        **dim_kwargs
            The keyword arguments form of ``dim``.
            One of dim or dim_kwargs must be provided.

        Returns
        -------
        stacked : Dataset
            Dataset with stacked data.

        See Also
        --------
        Dataset.unstack
        """
        dim = either_dict_or_kwargs(dim, dim_kwargs, "stack")
        result = self
        for new_dim, dims in dim.items():
            result = result._stack_once(dims, new_dim, index_cls, create_index)
        return result

    def to_stacked_array(
        self,
        new_dim: Hashable,
        sample_dims: Collection[Hashable],
        variable_dim: Hashable = "variable",
        name: Hashable | None = None,
    ) -> DataArray:
        """Combine variables of differing dimensionality into a DataArray
        without broadcasting.

        This method is similar to Dataset.to_dataarray but does not broadcast the
        variables.

        Parameters
        ----------
        new_dim : hashable
            Name of the new stacked coordinate
        sample_dims : Collection of hashables
            List of dimensions that **will not** be stacked. Each array in the
            dataset must share these dimensions. For machine learning
            applications, these define the dimensions over which samples are
            drawn.
        variable_dim : hashable, default: "variable"
            Name of the level in the stacked coordinate which corresponds to
            the variables.
        name : hashable, optional
            Name of the new data array.

        Returns
        -------
        stacked : DataArray
            DataArray with the specified dimensions and data variables
            stacked together. The stacked coordinate is named ``new_dim``
            and represented by a MultiIndex object with a level containing the
            data variable names. The name of this level is controlled using
            the ``variable_dim`` argument.

        See Also
        --------
        Dataset.to_dataarray
        Dataset.stack
        DataArray.to_unstacked_dataset

        Examples
        --------
        >>> data = xr.Dataset(
        ...     data_vars={
        ...         "a": (("x", "y"), [[0, 1, 2], [3, 4, 5]]),
        ...         "b": ("x", [6, 7]),
        ...     },
        ...     coords={"y": ["u", "v", "w"]},
        ... )

        >>> data
        <xarray.Dataset> Size: 76B
        Dimensions:  (x: 2, y: 3)
        Coordinates:
          * y        (y) <U1 12B 'u' 'v' 'w'
        Dimensions without coordinates: x
        Data variables:
            a        (x, y) int64 48B 0 1 2 3 4 5
            b        (x) int64 16B 6 7

        >>> data.to_stacked_array("z", sample_dims=["x"])
        <xarray.DataArray 'a' (x: 2, z: 4)> Size: 64B
        array([[0, 1, 2, 6],
               [3, 4, 5, 7]])
        Coordinates:
          * z         (z) object 32B MultiIndex
          * variable  (z) <U1 16B 'a' 'a' 'a' 'b'
          * y         (z) object 32B 'u' 'v' 'w' nan
        Dimensions without coordinates: x

        """
        from xarray.structure.concat import concat

        # add stacking dims by order of appearance
        stacking_dims_list: list[Hashable] = []
        for da in self.data_vars.values():
            for dim in da.dims:
                if dim not in sample_dims and dim not in stacking_dims_list:
                    stacking_dims_list.append(dim)
        stacking_dims = tuple(stacking_dims_list)

        for key, da in self.data_vars.items():
            missing_sample_dims = set(sample_dims) - set(da.dims)
            if missing_sample_dims:
                raise ValueError(
                    "Variables in the dataset must contain all ``sample_dims`` "
                    f"({sample_dims!r}) but '{key}' misses {sorted(map(str, missing_sample_dims))}"
                )

        def stack_dataarray(da):
            # add missing dims/ coords and the name of the variable

            missing_stack_coords = {variable_dim: da.name}
            for dim in set(stacking_dims) - set(da.dims):
                missing_stack_coords[dim] = None

            missing_stack_dims = list(missing_stack_coords)

            return (
                da.assign_coords(**missing_stack_coords)
                .expand_dims(missing_stack_dims)
                .stack({new_dim: (variable_dim,) + stacking_dims})
            )

        # concatenate the arrays
        stackable_vars = [stack_dataarray(da) for da in self.data_vars.values()]
        data_array = concat(
            stackable_vars,
            dim=new_dim,
            data_vars="all",
            coords="different",
            compat="equals",
            join="outer",
        )

        if name is not None:
            data_array.name = name

        return data_array

    def _unstack_once(
        self,
        dim: Hashable,
        index_and_vars: tuple[Index, dict[Hashable, Variable]],
        fill_value,
        sparse: bool = False,
    ) -> Self:
        index, index_vars = index_and_vars
        variables: dict[Hashable, Variable] = {}
        indexes = {k: v for k, v in self._indexes.items() if k != dim}

        new_indexes, clean_index = index.unstack()
        indexes.update(new_indexes)

        for idx in new_indexes.values():
            variables.update(idx.create_variables(index_vars))

        for name, var in self.variables.items():
            if name not in index_vars:
                if dim in var.dims:
                    if isinstance(fill_value, Mapping):
                        fill_value_ = fill_value[name]
                    else:
                        fill_value_ = fill_value

                    variables[name] = var._unstack_once(
                        index=clean_index,
                        dim=dim,
                        fill_value=fill_value_,
                        sparse=sparse,
                    )
                else:
                    variables[name] = var

        coord_names = set(self._coord_names) - {dim} | set(new_indexes)

        return self._replace_with_new_dims(
            variables, coord_names=coord_names, indexes=indexes
        )

    def _unstack_full_reindex(
        self,
        dim: Hashable,
        index_and_vars: tuple[Index, dict[Hashable, Variable]],
        fill_value,
        sparse: bool,
    ) -> Self:
        index, index_vars = index_and_vars
        variables: dict[Hashable, Variable] = {}
        indexes = {k: v for k, v in self._indexes.items() if k != dim}

        new_indexes, clean_index = index.unstack()
        indexes.update(new_indexes)

        new_index_variables = {}
        for idx in new_indexes.values():
            new_index_variables.update(idx.create_variables(index_vars))

        new_dim_sizes = {k: v.size for k, v in new_index_variables.items()}
        variables.update(new_index_variables)

        # take a shortcut in case the MultiIndex was not modified.
        full_idx = pd.MultiIndex.from_product(
            clean_index.levels, names=clean_index.names
        )
        if clean_index.equals(full_idx):
            obj = self
        else:
            # TODO: we may depreciate implicit re-indexing with a pandas.MultiIndex
            xr_full_idx = PandasMultiIndex(full_idx, dim)
            indexers = Indexes(
                dict.fromkeys(index_vars, xr_full_idx),
                xr_full_idx.create_variables(index_vars),
            )
            obj = self._reindex(
                indexers, copy=False, fill_value=fill_value, sparse=sparse
            )

        for name, var in obj.variables.items():
            if name not in index_vars:
                if dim in var.dims:
                    variables[name] = var.unstack({dim: new_dim_sizes})
                else:
                    variables[name] = var

        coord_names = set(self._coord_names) - {dim} | set(new_dim_sizes)

        return self._replace_with_new_dims(
            variables, coord_names=coord_names, indexes=indexes
        )

    def unstack(
        self,
        dim: Dims = None,
        *,
        fill_value: Any = xrdtypes.NA,
        sparse: bool = False,
    ) -> Self:
        """
        Unstack existing dimensions corresponding to MultiIndexes into
        multiple new dimensions.

        New dimensions will be added at the end.

        Parameters
        ----------
        dim : str, Iterable of Hashable or None, optional
            Dimension(s) over which to unstack. By default unstacks all
            MultiIndexes.
        fill_value : scalar or dict-like, default: nan
            value to be filled. If a dict-like, maps variable names to
            fill values. If not provided or if the dict-like does not
            contain all variables, the dtype's NA value will be used.
        sparse : bool, default: False
            use sparse-array if True

        Returns
        -------
        unstacked : Dataset
            Dataset with unstacked data.

        See Also
        --------
        Dataset.stack
        """

        if dim is None:
            dims = list(self.dims)
        else:
            if isinstance(dim, str) or not isinstance(dim, Iterable):
                dims = [dim]
            else:
                dims = list(dim)

            missing_dims = set(dims) - set(self.dims)
            if missing_dims:
                raise ValueError(
                    f"Dimensions {tuple(missing_dims)} not found in data dimensions {tuple(self.dims)}"
                )

        # each specified dimension must have exactly one multi-index
        stacked_indexes: dict[Any, tuple[Index, dict[Hashable, Variable]]] = {}
        for d in dims:
            idx, idx_vars = self._get_stack_index(d, multi=True)
            if idx is not None:
                stacked_indexes[d] = idx, idx_vars

        if dim is None:
            dims = list(stacked_indexes)
        else:
            non_multi_dims = set(dims) - set(stacked_indexes)
            if non_multi_dims:
                raise ValueError(
                    "cannot unstack dimensions that do not "
                    f"have exactly one multi-index: {tuple(non_multi_dims)}"
                )

        result = self.copy(deep=False)

        # we want to avoid allocating an object-dtype ndarray for a MultiIndex,
        # so we can't just access self.variables[v].data for every variable.
        # We only check the non-index variables.
        # https://github.com/pydata/xarray/issues/5902
        nonindexes = [
            self.variables[k] for k in set(self.variables) - set(self._indexes)
        ]
        # Notes for each of these cases:
        # 1. Dask arrays don't support assignment by index, which the fast unstack
        #    function requires.
        #    https://github.com/pydata/xarray/pull/4746#issuecomment-753282125
        # 2. Sparse doesn't currently support (though we could special-case it)
        #    https://github.com/pydata/sparse/issues/422
        # 3. pint requires checking if it's a NumPy array until
        #    https://github.com/pydata/xarray/pull/4751 is resolved,
        #    Once that is resolved, explicitly exclude pint arrays.
        #    pint doesn't implement `np.full_like` in a way that's
        #    currently compatible.
        sparse_array_type = array_type("sparse")
        needs_full_reindex = any(
            is_duck_dask_array(v.data)
            or isinstance(v.data, sparse_array_type)
            or not isinstance(v.data, np.ndarray)
            for v in nonindexes
        )

        for d in dims:
            if needs_full_reindex:
                result = result._unstack_full_reindex(
                    d, stacked_indexes[d], fill_value, sparse
                )
            else:
                result = result._unstack_once(d, stacked_indexes[d], fill_value, sparse)
        return result

    def update(self, other: CoercibleMapping) -> Self:
        """Update this dataset's variables with those from another dataset.

        Just like :py:meth:`dict.update` this is a in-place operation.
        For a non-inplace version, see :py:meth:`Dataset.merge`.

        Parameters
        ----------
        other : Dataset or mapping
            Variables with which to update this dataset. One of:

            - Dataset
            - mapping {var name: DataArray}
            - mapping {var name: Variable}
            - mapping {var name: (dimension name, array-like)}
            - mapping {var name: (tuple of dimension names, array-like)}

        Returns
        -------
        updated : Dataset
            Updated dataset. Note that since the update is in-place this is the input
            dataset.

            It is deprecated since version 0.17 and scheduled to be removed in 0.21.

        Raises
        ------
        ValueError
            If any dimensions would have inconsistent sizes in the updated
            dataset.

        See Also
        --------
        Dataset.assign
        Dataset.merge
        """
        merge_result = dataset_update_method(self, other)
        return self._replace(inplace=True, **merge_result._asdict())

    def merge(
        self,
        other: CoercibleMapping | DataArray,
        overwrite_vars: Hashable | Iterable[Hashable] = frozenset(),
        compat: CompatOptions | CombineKwargDefault = _COMPAT_DEFAULT,
        join: JoinOptions | CombineKwargDefault = _JOIN_DEFAULT,
        fill_value: Any = xrdtypes.NA,
        combine_attrs: CombineAttrsOptions = "override",
    ) -> Self:
        """Merge the arrays of two datasets into a single dataset.

        This method generally does not allow for overriding data, with the
        exception of attributes, which are ignored on the second dataset.
        Variables with the same name are checked for conflicts via the equals
        or identical methods.

        Parameters
        ----------
        other : Dataset or mapping
            Dataset or variables to merge with this dataset.
        overwrite_vars : hashable or iterable of hashable, optional
            If provided, update variables of these name(s) without checking for
            conflicts in this dataset.
        compat : {"identical", "equals", "broadcast_equals", \
                  "no_conflicts", "override", "minimal"}, default: "no_conflicts"
            String indicating how to compare variables of the same name for
            potential conflicts:

            - 'identical': all values, dimensions and attributes must be the
              same.
            - 'equals': all values and dimensions must be the same.
            - 'broadcast_equals': all values must be equal when variables are
              broadcast against each other to ensure common dimensions.
            - 'no_conflicts': only values which are not null in both datasets
              must be equal. The returned dataset then contains the combination
              of all non-null values.
            - 'override': skip comparing and pick variable from first dataset
            - 'minimal': drop conflicting coordinates

        join : {"outer", "inner", "left", "right", "exact", "override"}, \
               default: "outer"
            Method for joining ``self`` and ``other`` along shared dimensions:

            - 'outer': use the union of the indexes
            - 'inner': use the intersection of the indexes
            - 'left': use indexes from ``self``
            - 'right': use indexes from ``other``
            - 'exact': error instead of aligning non-equal indexes
            - 'override': use indexes from ``self`` that are the same size
              as those of ``other`` in that dimension

        fill_value : scalar or dict-like, optional
            Value to use for newly missing values. If a dict-like, maps
            variable names (including coordinates) to fill values.
        combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts", \
                         "override"} or callable, default: "override"
            A callable or a string indicating how to combine attrs of the objects being
            merged:

            - "drop": empty attrs on returned Dataset.
            - "identical": all attrs must be the same on every object.
            - "no_conflicts": attrs from all objects are combined, any that have
              the same name must also have the same value.
            - "drop_conflicts": attrs from all objects are combined, any that have
              the same name but different values are dropped.
            - "override": skip comparing and copy attrs from the first dataset to
              the result.

            If a callable, it must expect a sequence of ``attrs`` dicts and a context object
            as its only parameters.

        Returns
        -------
        merged : Dataset
            Merged dataset.

        Raises
        ------
        MergeError
            If any variables conflict (see ``compat``).

        See Also
        --------
        Dataset.update
        """
        from xarray.core.dataarray import DataArray

        other = other.to_dataset() if isinstance(other, DataArray) else other
        merge_result = dataset_merge_method(
            self,
            other,
            overwrite_vars=overwrite_vars,
            compat=compat,
            join=join,
            fill_value=fill_value,
            combine_attrs=combine_attrs,
        )
        return self._replace(**merge_result._asdict())

    def _assert_all_in_dataset(
        self, names: Iterable[Hashable], virtual_okay: bool = False
    ) -> None:
        bad_names = set(names) - set(self._variables)
        if virtual_okay:
            bad_names -= self.virtual_variables
        if bad_names:
            ordered_bad_names = [name for name in names if name in bad_names]
            raise ValueError(
                f"These variables cannot be found in this dataset: {ordered_bad_names}"
            )

    def drop_vars(
        self,
        names: str | Iterable[Hashable] | Callable[[Self], str | Iterable[Hashable]],
        *,
        errors: ErrorOptions = "raise",
    ) -> Self:
        """Drop variables from this dataset.

        Parameters
        ----------
        names : Hashable or iterable of Hashable or Callable
            Name(s) of variables to drop. If a Callable, this object is passed as its
            only argument and its result is used.
        errors : {"raise", "ignore"}, default: "raise"
            If 'raise', raises a ValueError error if any of the variable
            passed are not in the dataset. If 'ignore', any given names that are in the
            dataset are dropped and no error is raised.

        Examples
        --------

        >>> dataset = xr.Dataset(
        ...     {
        ...         "temperature": (
        ...             ["time", "latitude", "longitude"],
        ...             [[[25.5, 26.3], [27.1, 28.0]]],
        ...         ),
        ...         "humidity": (
        ...             ["time", "latitude", "longitude"],
        ...             [[[65.0, 63.8], [58.2, 59.6]]],
        ...         ),
        ...         "wind_speed": (
        ...             ["time", "latitude", "longitude"],
        ...             [[[10.2, 8.5], [12.1, 9.8]]],
        ...         ),
        ...     },
        ...     coords={
        ...         "time": pd.date_range("2023-07-01", periods=1),
        ...         "latitude": [40.0, 40.2],
        ...         "longitude": [-75.0, -74.8],
        ...     },
        ... )
        >>> dataset
        <xarray.Dataset> Size: 136B
        Dimensions:      (time: 1, latitude: 2, longitude: 2)
        Coordinates:
          * time         (time) datetime64[ns] 8B 2023-07-01
          * latitude     (latitude) float64 16B 40.0 40.2
          * longitude    (longitude) float64 16B -75.0 -74.8
        Data variables:
            temperature  (time, latitude, longitude) float64 32B 25.5 26.3 27.1 28.0
            humidity     (time, latitude, longitude) float64 32B 65.0 63.8 58.2 59.6
            wind_speed   (time, latitude, longitude) float64 32B 10.2 8.5 12.1 9.8

        Drop the 'humidity' variable

        >>> dataset.drop_vars(["humidity"])
        <xarray.Dataset> Size: 104B
        Dimensions:      (time: 1, latitude: 2, longitude: 2)
        Coordinates:
          * time         (time) datetime64[ns] 8B 2023-07-01
          * latitude     (latitude) float64 16B 40.0 40.2
          * longitude    (longitude) float64 16B -75.0 -74.8
        Data variables:
            temperature  (time, latitude, longitude) float64 32B 25.5 26.3 27.1 28.0
            wind_speed   (time, latitude, longitude) float64 32B 10.2 8.5 12.1 9.8

        Drop the 'humidity', 'temperature' variables

        >>> dataset.drop_vars(["humidity", "temperature"])
        <xarray.Dataset> Size: 72B
        Dimensions:     (time: 1, latitude: 2, longitude: 2)
        Coordinates:
          * time        (time) datetime64[ns] 8B 2023-07-01
          * latitude    (latitude) float64 16B 40.0 40.2
          * longitude   (longitude) float64 16B -75.0 -74.8
        Data variables:
            wind_speed  (time, latitude, longitude) float64 32B 10.2 8.5 12.1 9.8

        Drop all indexes

        >>> dataset.drop_vars(lambda x: x.indexes)
        <xarray.Dataset> Size: 96B
        Dimensions:      (time: 1, latitude: 2, longitude: 2)
        Dimensions without coordinates: time, latitude, longitude
        Data variables:
            temperature  (time, latitude, longitude) float64 32B 25.5 26.3 27.1 28.0
            humidity     (time, latitude, longitude) float64 32B 65.0 63.8 58.2 59.6
            wind_speed   (time, latitude, longitude) float64 32B 10.2 8.5 12.1 9.8

        Attempt to drop non-existent variable with errors="ignore"

        >>> dataset.drop_vars(["pressure"], errors="ignore")
        <xarray.Dataset> Size: 136B
        Dimensions:      (time: 1, latitude: 2, longitude: 2)
        Coordinates:
          * time         (time) datetime64[ns] 8B 2023-07-01
          * latitude     (latitude) float64 16B 40.0 40.2
          * longitude    (longitude) float64 16B -75.0 -74.8
        Data variables:
            temperature  (time, latitude, longitude) float64 32B 25.5 26.3 27.1 28.0
            humidity     (time, latitude, longitude) float64 32B 65.0 63.8 58.2 59.6
            wind_speed   (time, latitude, longitude) float64 32B 10.2 8.5 12.1 9.8

        Attempt to drop non-existent variable with errors="raise"

        >>> dataset.drop_vars(["pressure"], errors="raise")
        Traceback (most recent call last):
        ValueError: These variables cannot be found in this dataset: ['pressure']

        Raises
        ------
        ValueError
             Raised if you attempt to drop a variable which is not present, and the kwarg ``errors='raise'``.

        Returns
        -------
        dropped : Dataset

        See Also
        --------
        DataArray.drop_vars

        """
        if callable(names):
            names = names(self)
        # the Iterable check is required for mypy
        if is_scalar(names) or not isinstance(names, Iterable):
            names_set = {names}
        else:
            names_set = set(names)
        if errors == "raise":
            self._assert_all_in_dataset(names_set)

        # GH6505
        other_names = set()
        for var in names_set:
            maybe_midx = self._indexes.get(var, None)
            if isinstance(maybe_midx, PandasMultiIndex):
                idx_coord_names = set(list(maybe_midx.index.names) + [maybe_midx.dim])
                idx_other_names = idx_coord_names - set(names_set)
                other_names.update(idx_other_names)
        if other_names:
            names_set |= set(other_names)
            emit_user_level_warning(
                f"Deleting a single level of a MultiIndex is deprecated. Previously, this deleted all levels of a MultiIndex. "
                f"Please also drop the following variables: {other_names!r} to avoid an error in the future.",
                DeprecationWarning,
            )

        assert_no_index_corrupted(self.xindexes, names_set)

        variables = {k: v for k, v in self._variables.items() if k not in names_set}
        coord_names = {k for k in self._coord_names if k in variables}
        indexes = {k: v for k, v in self._indexes.items() if k not in names_set}
        return self._replace_with_new_dims(
            variables, coord_names=coord_names, indexes=indexes
        )

    def drop_indexes(
        self,
        coord_names: Hashable | Iterable[Hashable],
        *,
        errors: ErrorOptions = "raise",
    ) -> Self:
        """Drop the indexes assigned to the given coordinates.

        Parameters
        ----------
        coord_names : hashable or iterable of hashable
            Name(s) of the coordinate(s) for which to drop the index.
        errors : {"raise", "ignore"}, default: "raise"
            If 'raise', raises a ValueError error if any of the coordinates
            passed have no index or are not in the dataset.
            If 'ignore', no error is raised.

        Returns
        -------
        dropped : Dataset
            A new dataset with dropped indexes.

        """
        # the Iterable check is required for mypy
        if is_scalar(coord_names) or not isinstance(coord_names, Iterable):
            coord_names = {coord_names}
        else:
            coord_names = set(coord_names)

        if errors == "raise":
            invalid_coords = coord_names - self._coord_names
            if invalid_coords:
                raise ValueError(
                    f"The coordinates {tuple(invalid_coords)} are not found in the "
                    f"dataset coordinates {tuple(self.coords.keys())}"
                )

            unindexed_coords = set(coord_names) - set(self._indexes)
            if unindexed_coords:
                raise ValueError(
                    f"those coordinates do not have an index: {unindexed_coords}"
                )

        assert_no_index_corrupted(self.xindexes, coord_names, action="remove index(es)")

        variables = {}
        for name, var in self._variables.items():
            if name in coord_names:
                variables[name] = var.to_base_variable()
            else:
                variables[name] = var

        indexes = {k: v for k, v in self._indexes.items() if k not in coord_names}

        return self._replace(variables=variables, indexes=indexes)

    def drop(
        self,
        labels=None,
        dim=None,
        *,
        errors: ErrorOptions = "raise",
        **labels_kwargs,
    ) -> Self:
        """Backward compatible method based on `drop_vars` and `drop_sel`

        Using either `drop_vars` or `drop_sel` is encouraged

        See Also
        --------
        Dataset.drop_vars
        Dataset.drop_sel
        """
        if errors not in ["raise", "ignore"]:
            raise ValueError('errors must be either "raise" or "ignore"')

        if is_dict_like(labels) and not isinstance(labels, dict):
            emit_user_level_warning(
                "dropping coordinates using `drop` is deprecated; use drop_vars.",
                DeprecationWarning,
            )
            return self.drop_vars(labels, errors=errors)

        if labels_kwargs or isinstance(labels, dict):
            if dim is not None:
                raise ValueError("cannot specify dim and dict-like arguments.")
            labels = either_dict_or_kwargs(labels, labels_kwargs, "drop")

        if dim is None and (is_scalar(labels) or isinstance(labels, Iterable)):
            emit_user_level_warning(
                "dropping variables using `drop` is deprecated; use drop_vars.",
                DeprecationWarning,
            )
            # for mypy
            if is_scalar(labels):
                labels = [labels]
            return self.drop_vars(labels, errors=errors)
        if dim is not None:
            warnings.warn(
                "dropping labels using list-like labels is deprecated; using "
                "dict-like arguments with `drop_sel`, e.g. `ds.drop_sel(dim=[labels]).",
                DeprecationWarning,
                stacklevel=2,
            )
            return self.drop_sel({dim: labels}, errors=errors, **labels_kwargs)

        emit_user_level_warning(
            "dropping labels using `drop` is deprecated; use `drop_sel` instead.",
            DeprecationWarning,
        )
        return self.drop_sel(labels, errors=errors)

    def drop_sel(
        self, labels=None, *, errors: ErrorOptions = "raise", **labels_kwargs
    ) -> Self:
        """Drop index labels from this dataset.

        Parameters
        ----------
        labels : mapping of hashable to Any
            Index labels to drop
        errors : {"raise", "ignore"}, default: "raise"
            If 'raise', raises a ValueError error if
            any of the index labels passed are not
            in the dataset. If 'ignore', any given labels that are in the
            dataset are dropped and no error is raised.
        **labels_kwargs : {dim: label, ...}, optional
            The keyword arguments form of ``dim`` and ``labels``

        Returns
        -------
        dropped : Dataset

        Examples
        --------
        >>> data = np.arange(6).reshape(2, 3)
        >>> labels = ["a", "b", "c"]
        >>> ds = xr.Dataset({"A": (["x", "y"], data), "y": labels})
        >>> ds
        <xarray.Dataset> Size: 60B
        Dimensions:  (x: 2, y: 3)
        Coordinates:
          * y        (y) <U1 12B 'a' 'b' 'c'
        Dimensions without coordinates: x
        Data variables:
            A        (x, y) int64 48B 0 1 2 3 4 5
        >>> ds.drop_sel(y=["a", "c"])
        <xarray.Dataset> Size: 20B
        Dimensions:  (x: 2, y: 1)
        Coordinates:
          * y        (y) <U1 4B 'b'
        Dimensions without coordinates: x
        Data variables:
            A        (x, y) int64 16B 1 4
        >>> ds.drop_sel(y="b")
        <xarray.Dataset> Size: 40B
        Dimensions:  (x: 2, y: 2)
        Coordinates:
          * y        (y) <U1 8B 'a' 'c'
        Dimensions without coordinates: x
        Data variables:
            A        (x, y) int64 32B 0 2 3 5
        """
        if errors not in ["raise", "ignore"]:
            raise ValueError('errors must be either "raise" or "ignore"')

        labels = either_dict_or_kwargs(labels, labels_kwargs, "drop_sel")

        ds = self
        for dim, labels_for_dim in labels.items():
            # Don't cast to set, as it would harm performance when labels
            # is a large numpy array
            if utils.is_scalar(labels_for_dim):
                labels_for_dim = [labels_for_dim]
            labels_for_dim = np.asarray(labels_for_dim)
            try:
                index = self.get_index(dim)
            except KeyError as err:
                raise ValueError(
                    f"dimension {dim!r} does not have coordinate labels"
                ) from err
            new_index = index.drop(labels_for_dim, errors=errors)
            ds = ds.loc[{dim: new_index}]
        return ds

    def drop_isel(self, indexers=None, **indexers_kwargs) -> Self:
        """Drop index positions from this Dataset.

        Parameters
        ----------
        indexers : mapping of hashable to Any
            Index locations to drop
        **indexers_kwargs : {dim: position, ...}, optional
            The keyword arguments form of ``dim`` and ``positions``

        Returns
        -------
        dropped : Dataset

        Raises
        ------
        IndexError

        Examples
        --------
        >>> data = np.arange(6).reshape(2, 3)
        >>> labels = ["a", "b", "c"]
        >>> ds = xr.Dataset({"A": (["x", "y"], data), "y": labels})
        >>> ds
        <xarray.Dataset> Size: 60B
        Dimensions:  (x: 2, y: 3)
        Coordinates:
          * y        (y) <U1 12B 'a' 'b' 'c'
        Dimensions without coordinates: x
        Data variables:
            A        (x, y) int64 48B 0 1 2 3 4 5
        >>> ds.drop_isel(y=[0, 2])
        <xarray.Dataset> Size: 20B
        Dimensions:  (x: 2, y: 1)
        Coordinates:
          * y        (y) <U1 4B 'b'
        Dimensions without coordinates: x
        Data variables:
            A        (x, y) int64 16B 1 4
        >>> ds.drop_isel(y=1)
        <xarray.Dataset> Size: 40B
        Dimensions:  (x: 2, y: 2)
        Coordinates:
          * y        (y) <U1 8B 'a' 'c'
        Dimensions without coordinates: x
        Data variables:
            A        (x, y) int64 32B 0 2 3 5
        """

        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "drop_isel")

        ds = self
        dimension_index = {}
        for dim, pos_for_dim in indexers.items():
            # Don't cast to set, as it would harm performance when labels
            # is a large numpy array
            if utils.is_scalar(pos_for_dim):
                pos_for_dim = [pos_for_dim]
            pos_for_dim = np.asarray(pos_for_dim)
            index = self.get_index(dim)
            new_index = index.delete(pos_for_dim)
            dimension_index[dim] = new_index
        ds = ds.loc[dimension_index]
        return ds

    def drop_dims(
        self,
        drop_dims: str | Iterable[Hashable],
        *,
        errors: ErrorOptions = "raise",
    ) -> Self:
        """Drop dimensions and associated variables from this dataset.

        Parameters
        ----------
        drop_dims : str or Iterable of Hashable
            Dimension or dimensions to drop.
        errors : {"raise", "ignore"}, default: "raise"
            If 'raise', raises a ValueError error if any of the
            dimensions passed are not in the dataset. If 'ignore', any given
            dimensions that are in the dataset are dropped and no error is raised.

        Returns
        -------
        obj : Dataset
            The dataset without the given dimensions (or any variables
            containing those dimensions).
        """
        if errors not in ["raise", "ignore"]:
            raise ValueError('errors must be either "raise" or "ignore"')

        if isinstance(drop_dims, str) or not isinstance(drop_dims, Iterable):
            drop_dims = {drop_dims}
        else:
            drop_dims = set(drop_dims)

        if errors == "raise":
            missing_dims = drop_dims - set(self.dims)
            if missing_dims:
                raise ValueError(
                    f"Dimensions {tuple(missing_dims)} not found in data dimensions {tuple(self.dims)}"
                )

        drop_vars = {k for k, v in self._variables.items() if set(v.dims) & drop_dims}
        return self.drop_vars(drop_vars)

    @deprecate_dims
    def transpose(
        self,
        *dim: Hashable,
        missing_dims: ErrorOptionsWithWarn = "raise",
    ) -> Self:
        """Return a new Dataset object with all array dimensions transposed.

        Although the order of dimensions on each array will change, the dataset
        dimensions themselves will remain in fixed (sorted) order.

        Parameters
        ----------
        *dim : hashable, optional
            By default, reverse the dimensions on each array. Otherwise,
            reorder the dimensions to this order.
        missing_dims : {"raise", "warn", "ignore"}, default: "raise"
            What to do if dimensions that should be selected from are not present in the
            Dataset:
            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        Returns
        -------
        transposed : Dataset
            Each array in the dataset (including) coordinates will be
            transposed to the given order.

        Notes
        -----
        This operation returns a view of each array's data. It is
        lazy for dask-backed DataArrays but not for numpy-backed DataArrays
        -- the data will be fully loaded into memory.

        See Also
        --------
        numpy.transpose
        DataArray.transpose
        """
        # Raise error if list is passed as dim
        if (len(dim) > 0) and (isinstance(dim[0], list)):
            list_fix = [f"{x!r}" if isinstance(x, str) else f"{x}" for x in dim[0]]
            raise TypeError(
                f"transpose requires dim to be passed as multiple arguments. Expected `{', '.join(list_fix)}`. Received `{dim[0]}` instead"
            )

        # Use infix_dims to check once for missing dimensions
        if len(dim) != 0:
            _ = list(infix_dims(dim, self.dims, missing_dims))

        ds = self.copy()
        for name, var in self._variables.items():
            var_dims = tuple(d for d in dim if d in (var.dims + (...,)))
            ds._variables[name] = var.transpose(*var_dims)
        return ds

    def dropna(
        self,
        dim: Hashable,
        *,
        how: Literal["any", "all"] = "any",
        thresh: int | None = None,
        subset: Iterable[Hashable] | None = None,
    ) -> Self:
        """Returns a new dataset with dropped labels for missing values along
        the provided dimension.

        Parameters
        ----------
        dim : hashable
            Dimension along which to drop missing values. Dropping along
            multiple dimensions simultaneously is not yet supported.
        how : {"any", "all"}, default: "any"
            - any : if any NA values are present, drop that label
            - all : if all values are NA, drop that label

        thresh : int or None, optional
            If supplied, require this many non-NA values (summed over all the subset variables).
        subset : iterable of hashable or None, optional
            Which variables to check for missing values. By default, all
            variables in the dataset are checked.

        Examples
        --------
        >>> dataset = xr.Dataset(
        ...     {
        ...         "temperature": (
        ...             ["time", "location"],
        ...             [[23.4, 24.1], [np.nan, 22.1], [21.8, 24.2], [20.5, 25.3]],
        ...         )
        ...     },
        ...     coords={"time": [1, 2, 3, 4], "location": ["A", "B"]},
        ... )
        >>> dataset
        <xarray.Dataset> Size: 104B
        Dimensions:      (time: 4, location: 2)
        Coordinates:
          * time         (time) int64 32B 1 2 3 4
          * location     (location) <U1 8B 'A' 'B'
        Data variables:
            temperature  (time, location) float64 64B 23.4 24.1 nan ... 24.2 20.5 25.3

        Drop NaN values from the dataset

        >>> dataset.dropna(dim="time")
        <xarray.Dataset> Size: 80B
        Dimensions:      (time: 3, location: 2)
        Coordinates:
          * time         (time) int64 24B 1 3 4
          * location     (location) <U1 8B 'A' 'B'
        Data variables:
            temperature  (time, location) float64 48B 23.4 24.1 21.8 24.2 20.5 25.3

        Drop labels with any NaN values

        >>> dataset.dropna(dim="time", how="any")
        <xarray.Dataset> Size: 80B
        Dimensions:      (time: 3, location: 2)
        Coordinates:
          * time         (time) int64 24B 1 3 4
          * location     (location) <U1 8B 'A' 'B'
        Data variables:
            temperature  (time, location) float64 48B 23.4 24.1 21.8 24.2 20.5 25.3

        Drop labels with all NAN values

        >>> dataset.dropna(dim="time", how="all")
        <xarray.Dataset> Size: 104B
        Dimensions:      (time: 4, location: 2)
        Coordinates:
          * time         (time) int64 32B 1 2 3 4
          * location     (location) <U1 8B 'A' 'B'
        Data variables:
            temperature  (time, location) float64 64B 23.4 24.1 nan ... 24.2 20.5 25.3

        Drop labels with less than 2 non-NA values

        >>> dataset.dropna(dim="time", thresh=2)
        <xarray.Dataset> Size: 80B
        Dimensions:      (time: 3, location: 2)
        Coordinates:
          * time         (time) int64 24B 1 3 4
          * location     (location) <U1 8B 'A' 'B'
        Data variables:
            temperature  (time, location) float64 48B 23.4 24.1 21.8 24.2 20.5 25.3

        Returns
        -------
        Dataset

        See Also
        --------
        DataArray.dropna
        pandas.DataFrame.dropna
        """
        # TODO: consider supporting multiple dimensions? Or not, given that
        # there are some ugly edge cases, e.g., pandas's dropna differs
        # depending on the order of the supplied axes.

        if dim not in self.dims:
            raise ValueError(
                f"Dimension {dim!r} not found in data dimensions {tuple(self.dims)}"
            )

        if subset is None:
            subset = iter(self.data_vars)

        count = np.zeros(self.sizes[dim], dtype=np.int64)
        size = np.int_(0)  # for type checking

        for k in subset:
            array = self._variables[k]
            if dim in array.dims:
                dims = [d for d in array.dims if d != dim]
                count += to_numpy(array.count(dims).data)
                size += math.prod([self.sizes[d] for d in dims])

        if thresh is not None:
            mask = count >= thresh
        elif how == "any":
            mask = count == size
        elif how == "all":
            mask = count > 0
        elif how is not None:
            raise ValueError(f"invalid how option: {how}")
        else:
            raise TypeError("must specify how or thresh")

        return self.isel({dim: mask})

    def fillna(self, value: Any) -> Self:
        """Fill missing values in this object.

        This operation follows the normal broadcasting and alignment rules that
        xarray uses for binary arithmetic, except the result is aligned to this
        object (``join='left'``) instead of aligned to the intersection of
        index coordinates (``join='inner'``).

        Parameters
        ----------
        value : scalar, ndarray, DataArray, dict or Dataset
            Used to fill all matching missing values in this dataset's data
            variables. Scalars, ndarrays or DataArrays arguments are used to
            fill all data with aligned coordinates (for DataArrays).
            Dictionaries or datasets match data variables and then align
            coordinates if necessary.

        Returns
        -------
        Dataset

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     {
        ...         "A": ("x", [np.nan, 2, np.nan, 0]),
        ...         "B": ("x", [3, 4, np.nan, 1]),
        ...         "C": ("x", [np.nan, np.nan, np.nan, 5]),
        ...         "D": ("x", [np.nan, 3, np.nan, 4]),
        ...     },
        ...     coords={"x": [0, 1, 2, 3]},
        ... )
        >>> ds
        <xarray.Dataset> Size: 160B
        Dimensions:  (x: 4)
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
        Data variables:
            A        (x) float64 32B nan 2.0 nan 0.0
            B        (x) float64 32B 3.0 4.0 nan 1.0
            C        (x) float64 32B nan nan nan 5.0
            D        (x) float64 32B nan 3.0 nan 4.0

        Replace all `NaN` values with 0s.

        >>> ds.fillna(0)
        <xarray.Dataset> Size: 160B
        Dimensions:  (x: 4)
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
        Data variables:
            A        (x) float64 32B 0.0 2.0 0.0 0.0
            B        (x) float64 32B 3.0 4.0 0.0 1.0
            C        (x) float64 32B 0.0 0.0 0.0 5.0
            D        (x) float64 32B 0.0 3.0 0.0 4.0

        Replace all `NaN` elements in column ‘A’, ‘B’, ‘C’, and ‘D’, with 0, 1, 2, and 3 respectively.

        >>> values = {"A": 0, "B": 1, "C": 2, "D": 3}
        >>> ds.fillna(value=values)
        <xarray.Dataset> Size: 160B
        Dimensions:  (x: 4)
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
        Data variables:
            A        (x) float64 32B 0.0 2.0 0.0 0.0
            B        (x) float64 32B 3.0 4.0 1.0 1.0
            C        (x) float64 32B 2.0 2.0 2.0 5.0
            D        (x) float64 32B 3.0 3.0 3.0 4.0
        """
        if utils.is_dict_like(value):
            value_keys = getattr(value, "data_vars", value).keys()
            if not set(value_keys) <= set(self.data_vars.keys()):
                raise ValueError(
                    "all variables in the argument to `fillna` "
                    "must be contained in the original dataset"
                )
        out = ops.fillna(self, value)
        return out

    def interpolate_na(
        self,
        dim: Hashable | None = None,
        method: InterpOptions = "linear",
        limit: int | None = None,
        use_coordinate: bool | Hashable = True,
        max_gap: (
            int
            | float
            | str
            | pd.Timedelta
            | np.timedelta64
            | datetime.timedelta
            | None
        ) = None,
        **kwargs: Any,
    ) -> Self:
        """Fill in NaNs by interpolating according to different methods.

        Parameters
        ----------
        dim : Hashable or None, optional
            Specifies the dimension along which to interpolate.
        method : {"linear", "nearest", "zero", "slinear", "quadratic", "cubic", "polynomial", \
            "barycentric", "krogh", "pchip", "spline", "akima"}, default: "linear"
            String indicating which method to use for interpolation:

            - 'linear': linear interpolation. Additional keyword
              arguments are passed to :py:func:`numpy.interp`
            - 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'polynomial':
              are passed to :py:func:`scipy.interpolate.interp1d`. If
              ``method='polynomial'``, the ``order`` keyword argument must also be
              provided.
            - 'barycentric', 'krogh', 'pchip', 'spline', 'akima': use their
              respective :py:class:`scipy.interpolate` classes.

        use_coordinate : bool or Hashable, default: True
            Specifies which index to use as the x values in the interpolation
            formulated as `y = f(x)`. If False, values are treated as if
            equally-spaced along ``dim``. If True, the IndexVariable `dim` is
            used. If ``use_coordinate`` is a string, it specifies the name of a
            coordinate variable to use as the index.
        limit : int, default: None
            Maximum number of consecutive NaNs to fill. Must be greater than 0
            or None for no limit. This filling is done regardless of the size of
            the gap in the data. To only interpolate over gaps less than a given length,
            see ``max_gap``.
        max_gap : int, float, str, pandas.Timedelta, numpy.timedelta64, datetime.timedelta \
            or None, default: None
            Maximum size of gap, a continuous sequence of NaNs, that will be filled.
            Use None for no limit. When interpolating along a datetime64 dimension
            and ``use_coordinate=True``, ``max_gap`` can be one of the following:

            - a string that is valid input for pandas.to_timedelta
            - a :py:class:`numpy.timedelta64` object
            - a :py:class:`pandas.Timedelta` object
            - a :py:class:`datetime.timedelta` object

            Otherwise, ``max_gap`` must be an int or a float. Use of ``max_gap`` with unlabeled
            dimensions has not been implemented yet. Gap length is defined as the difference
            between coordinate values at the first data point after a gap and the last value
            before a gap. For gaps at the beginning (end), gap length is defined as the difference
            between coordinate values at the first (last) valid data point and the first (last) NaN.
            For example, consider::

                <xarray.DataArray (x: 9)>
                array([nan, nan, nan,  1., nan, nan,  4., nan, nan])
                Coordinates:
                  * x        (x) int64 0 1 2 3 4 5 6 7 8

            The gap lengths are 3-0 = 3; 6-3 = 3; and 8-6 = 2 respectively
        **kwargs : dict, optional
            parameters passed verbatim to the underlying interpolation function

        Returns
        -------
        interpolated: Dataset
            Filled in Dataset.

        Warning
        --------
        When passing fill_value as a keyword argument with method="linear", it does not use
        ``numpy.interp`` but it uses ``scipy.interpolate.interp1d``, which provides the fill_value parameter.

        See Also
        --------
        numpy.interp
        scipy.interpolate

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     {
        ...         "A": ("x", [np.nan, 2, 3, np.nan, 0]),
        ...         "B": ("x", [3, 4, np.nan, 1, 7]),
        ...         "C": ("x", [np.nan, np.nan, np.nan, 5, 0]),
        ...         "D": ("x", [np.nan, 3, np.nan, -1, 4]),
        ...     },
        ...     coords={"x": [0, 1, 2, 3, 4]},
        ... )
        >>> ds
        <xarray.Dataset> Size: 200B
        Dimensions:  (x: 5)
        Coordinates:
          * x        (x) int64 40B 0 1 2 3 4
        Data variables:
            A        (x) float64 40B nan 2.0 3.0 nan 0.0
            B        (x) float64 40B 3.0 4.0 nan 1.0 7.0
            C        (x) float64 40B nan nan nan 5.0 0.0
            D        (x) float64 40B nan 3.0 nan -1.0 4.0

        >>> ds.interpolate_na(dim="x", method="linear")
        <xarray.Dataset> Size: 200B
        Dimensions:  (x: 5)
        Coordinates:
          * x        (x) int64 40B 0 1 2 3 4
        Data variables:
            A        (x) float64 40B nan 2.0 3.0 1.5 0.0
            B        (x) float64 40B 3.0 4.0 2.5 1.0 7.0
            C        (x) float64 40B nan nan nan 5.0 0.0
            D        (x) float64 40B nan 3.0 1.0 -1.0 4.0

        >>> ds.interpolate_na(dim="x", method="linear", fill_value="extrapolate")
        <xarray.Dataset> Size: 200B
        Dimensions:  (x: 5)
        Coordinates:
          * x        (x) int64 40B 0 1 2 3 4
        Data variables:
            A        (x) float64 40B 1.0 2.0 3.0 1.5 0.0
            B        (x) float64 40B 3.0 4.0 2.5 1.0 7.0
            C        (x) float64 40B 20.0 15.0 10.0 5.0 0.0
            D        (x) float64 40B 5.0 3.0 1.0 -1.0 4.0
        """
        from xarray.core.missing import _apply_over_vars_with_dim, interp_na

        new = _apply_over_vars_with_dim(
            interp_na,
            self,
            dim=dim,
            method=method,
            limit=limit,
            use_coordinate=use_coordinate,
            max_gap=max_gap,
            **kwargs,
        )
        return new

    def ffill(self, dim: Hashable, limit: int | None = None) -> Self:
        """Fill NaN values by propagating values forward

        *Requires bottleneck.*

        Parameters
        ----------
        dim : Hashable
            Specifies the dimension along which to propagate values when filling.
        limit : int or None, optional
            The maximum number of consecutive NaN values to forward fill. In
            other words, if there is a gap with more than this number of
            consecutive NaNs, it will only be partially filled. Must be greater
            than 0 or None for no limit. Must be None or greater than or equal
            to axis length if filling along chunked axes (dimensions).

        Examples
        --------
        >>> time = pd.date_range("2023-01-01", periods=10, freq="D")
        >>> data = np.array(
        ...     [1, np.nan, np.nan, np.nan, 5, np.nan, np.nan, 8, np.nan, 10]
        ... )
        >>> dataset = xr.Dataset({"data": (("time",), data)}, coords={"time": time})
        >>> dataset
        <xarray.Dataset> Size: 160B
        Dimensions:  (time: 10)
        Coordinates:
          * time     (time) datetime64[ns] 80B 2023-01-01 2023-01-02 ... 2023-01-10
        Data variables:
            data     (time) float64 80B 1.0 nan nan nan 5.0 nan nan 8.0 nan 10.0

        # Perform forward fill (ffill) on the dataset

        >>> dataset.ffill(dim="time")
        <xarray.Dataset> Size: 160B
        Dimensions:  (time: 10)
        Coordinates:
          * time     (time) datetime64[ns] 80B 2023-01-01 2023-01-02 ... 2023-01-10
        Data variables:
            data     (time) float64 80B 1.0 1.0 1.0 1.0 5.0 5.0 5.0 8.0 8.0 10.0

        # Limit the forward filling to a maximum of 2 consecutive NaN values

        >>> dataset.ffill(dim="time", limit=2)
        <xarray.Dataset> Size: 160B
        Dimensions:  (time: 10)
        Coordinates:
          * time     (time) datetime64[ns] 80B 2023-01-01 2023-01-02 ... 2023-01-10
        Data variables:
            data     (time) float64 80B 1.0 1.0 1.0 nan 5.0 5.0 5.0 8.0 8.0 10.0

        Returns
        -------
        Dataset

        See Also
        --------
        Dataset.bfill
        """
        from xarray.core.missing import _apply_over_vars_with_dim, ffill

        new = _apply_over_vars_with_dim(ffill, self, dim=dim, limit=limit)
        return new

    def bfill(self, dim: Hashable, limit: int | None = None) -> Self:
        """Fill NaN values by propagating values backward

        *Requires bottleneck.*

        Parameters
        ----------
        dim : Hashable
            Specifies the dimension along which to propagate values when
            filling.
        limit : int or None, optional
            The maximum number of consecutive NaN values to backward fill. In
            other words, if there is a gap with more than this number of
            consecutive NaNs, it will only be partially filled. Must be greater
            than 0 or None for no limit. Must be None or greater than or equal
            to axis length if filling along chunked axes (dimensions).

        Examples
        --------
        >>> time = pd.date_range("2023-01-01", periods=10, freq="D")
        >>> data = np.array(
        ...     [1, np.nan, np.nan, np.nan, 5, np.nan, np.nan, 8, np.nan, 10]
        ... )
        >>> dataset = xr.Dataset({"data": (("time",), data)}, coords={"time": time})
        >>> dataset
        <xarray.Dataset> Size: 160B
        Dimensions:  (time: 10)
        Coordinates:
          * time     (time) datetime64[ns] 80B 2023-01-01 2023-01-02 ... 2023-01-10
        Data variables:
            data     (time) float64 80B 1.0 nan nan nan 5.0 nan nan 8.0 nan 10.0

        # filled dataset, fills NaN values by propagating values backward

        >>> dataset.bfill(dim="time")
        <xarray.Dataset> Size: 160B
        Dimensions:  (time: 10)
        Coordinates:
          * time     (time) datetime64[ns] 80B 2023-01-01 2023-01-02 ... 2023-01-10
        Data variables:
            data     (time) float64 80B 1.0 5.0 5.0 5.0 5.0 8.0 8.0 8.0 10.0 10.0

        # Limit the backward filling to a maximum of 2 consecutive NaN values

        >>> dataset.bfill(dim="time", limit=2)
        <xarray.Dataset> Size: 160B
        Dimensions:  (time: 10)
        Coordinates:
          * time     (time) datetime64[ns] 80B 2023-01-01 2023-01-02 ... 2023-01-10
        Data variables:
            data     (time) float64 80B 1.0 nan 5.0 5.0 5.0 8.0 8.0 8.0 10.0 10.0

        Returns
        -------
        Dataset

        See Also
        --------
        Dataset.ffill
        """
        from xarray.core.missing import _apply_over_vars_with_dim, bfill

        new = _apply_over_vars_with_dim(bfill, self, dim=dim, limit=limit)
        return new

    def combine_first(self, other: Self) -> Self:
        """Combine two Datasets, default to data_vars of self.

        The new coordinates follow the normal broadcasting and alignment rules
        of ``join='outer'``.  Vacant cells in the expanded coordinates are
        filled with np.nan.

        Parameters
        ----------
        other : Dataset
            Used to fill all matching missing values in this array.

        Returns
        -------
        Dataset
        """
        out = ops.fillna(self, other, join="outer", dataset_join="outer")
        return out

    def reduce(
        self,
        func: Callable,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        numeric_only: bool = False,
        **kwargs: Any,
    ) -> Self:
        """Reduce this dataset by applying `func` along some dimension(s).

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `f(x, axis=axis, **kwargs)` to return the result of reducing an
            np.ndarray over an integer valued axis.
        dim : str, Iterable of Hashable or None, optional
            Dimension(s) over which to apply `func`. By default `func` is
            applied over all dimensions.
        keep_attrs : bool or None, optional
            If True, the dataset's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        keepdims : bool, default: False
            If True, the dimensions which are reduced are left in the result
            as dimensions of size one. Coordinates that use these dimensions
            are removed.
        numeric_only : bool, default: False
            If True, only apply ``func`` to variables with a numeric dtype.
        **kwargs : Any
            Additional keyword arguments passed on to ``func``.

        Returns
        -------
        reduced : Dataset
            Dataset with this object's DataArrays replaced with new DataArrays
            of summarized data and the indicated dimension(s) removed.

        Examples
        --------

        >>> dataset = xr.Dataset(
        ...     {
        ...         "math_scores": (
        ...             ["student", "test"],
        ...             [[90, 85, 92], [78, 80, 85], [95, 92, 98]],
        ...         ),
        ...         "english_scores": (
        ...             ["student", "test"],
        ...             [[88, 90, 92], [75, 82, 79], [93, 96, 91]],
        ...         ),
        ...     },
        ...     coords={
        ...         "student": ["Alice", "Bob", "Charlie"],
        ...         "test": ["Test 1", "Test 2", "Test 3"],
        ...     },
        ... )

        # Calculate the 75th percentile of math scores for each student using np.percentile

        >>> percentile_scores = dataset.reduce(np.percentile, q=75, dim="test")
        >>> percentile_scores
        <xarray.Dataset> Size: 132B
        Dimensions:         (student: 3)
        Coordinates:
          * student         (student) <U7 84B 'Alice' 'Bob' 'Charlie'
        Data variables:
            math_scores     (student) float64 24B 91.0 82.5 96.5
            english_scores  (student) float64 24B 91.0 80.5 94.5
        """
        if kwargs.get("axis") is not None:
            raise ValueError(
                "passing 'axis' to Dataset reduce methods is ambiguous."
                " Please use 'dim' instead."
            )

        dims = parse_dims_as_set(dim, set(self._dims.keys()))

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)

        variables: dict[Hashable, Variable] = {}
        for name, var in self._variables.items():
            reduce_dims = [d for d in var.dims if d in dims]
            if name in self.coords:
                if not reduce_dims:
                    variables[name] = var
            elif (
                # Some reduction functions (e.g. std, var) need to run on variables
                # that don't have the reduce dims: PR5393
                not pd.api.types.is_extension_array_dtype(var.dtype)  # noqa: TID251
                and (
                    not reduce_dims
                    or not numeric_only
                    or np.issubdtype(var.dtype, np.number)
                    or (var.dtype == np.bool_)
                )
            ):
                # prefer to aggregate over axis=None rather than
                # axis=(0, 1) if they will be equivalent, because
                # the former is often more efficient
                # keep single-element dims as list, to support Hashables
                reduce_maybe_single = (
                    None
                    if len(reduce_dims) == var.ndim and var.ndim != 1
                    else reduce_dims
                )
                variables[name] = var.reduce(
                    func,
                    dim=reduce_maybe_single,
                    keep_attrs=keep_attrs,
                    keepdims=keepdims,
                    **kwargs,
                )

        coord_names = {k for k in self.coords if k in variables}
        indexes = {k: v for k, v in self._indexes.items() if k in variables}
        attrs = self.attrs if keep_attrs else None
        return self._replace_with_new_dims(
            variables, coord_names=coord_names, attrs=attrs, indexes=indexes
        )

    def map(
        self,
        func: Callable,
        keep_attrs: bool | None = None,
        args: Iterable[Any] = (),
        **kwargs: Any,
    ) -> Self:
        """Apply a function to each data variable in this dataset

        Parameters
        ----------
        func : callable
            Function which can be called in the form `func(x, *args, **kwargs)`
            to transform each DataArray `x` in this dataset into another
            DataArray.
        keep_attrs : bool or None, optional
            If True, both the dataset's and variables' attributes (`attrs`) will be
            copied from the original objects to the new ones. If False, the new dataset
            and variables will be returned without copying the attributes.
        args : iterable, optional
            Positional arguments passed on to `func`.
        **kwargs : Any
            Keyword arguments passed on to `func`.

        Returns
        -------
        applied : Dataset
            Resulting dataset from applying ``func`` to each data variable.

        Examples
        --------
        >>> da = xr.DataArray(np.random.randn(2, 3))
        >>> ds = xr.Dataset({"foo": da, "bar": ("x", [-1, 2])})
        >>> ds
        <xarray.Dataset> Size: 64B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Dimensions without coordinates: dim_0, dim_1, x
        Data variables:
            foo      (dim_0, dim_1) float64 48B 1.764 0.4002 0.9787 2.241 1.868 -0.9773
            bar      (x) int64 16B -1 2
        >>> ds.map(np.fabs)
        <xarray.Dataset> Size: 64B
        Dimensions:  (dim_0: 2, dim_1: 3, x: 2)
        Dimensions without coordinates: dim_0, dim_1, x
        Data variables:
            foo      (dim_0, dim_1) float64 48B 1.764 0.4002 0.9787 2.241 1.868 0.9773
            bar      (x) float64 16B 1.0 2.0
        """
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)
        variables = {
            k: maybe_wrap_array(v, func(v, *args, **kwargs))
            for k, v in self.data_vars.items()
        }
        if keep_attrs:
            for k, v in variables.items():
                v._copy_attrs_from(self.data_vars[k])
        attrs = self.attrs if keep_attrs else None
        return type(self)(variables, attrs=attrs)

    def apply(
        self,
        func: Callable,
        keep_attrs: bool | None = None,
        args: Iterable[Any] = (),
        **kwargs: Any,
    ) -> Self:
        """
        Backward compatible implementation of ``map``

        See Also
        --------
        Dataset.map
        """
        warnings.warn(
            "Dataset.apply may be deprecated in the future. Using Dataset.map is encouraged",
            PendingDeprecationWarning,
            stacklevel=2,
        )
        return self.map(func, keep_attrs, args, **kwargs)

    def assign(
        self,
        variables: Mapping[Any, Any] | None = None,
        **variables_kwargs: Any,
    ) -> Self:
        """Assign new data variables to a Dataset, returning a new object
        with all the original variables in addition to the new ones.

        Parameters
        ----------
        variables : mapping of hashable to Any
            Mapping from variables names to the new values. If the new values
            are callable, they are computed on the Dataset and assigned to new
            data variables. If the values are not callable, (e.g. a DataArray,
            scalar, or array), they are simply assigned.
        **variables_kwargs
            The keyword arguments form of ``variables``.
            One of variables or variables_kwargs must be provided.

        Returns
        -------
        ds : Dataset
            A new Dataset with the new variables in addition to all the
            existing variables.

        Notes
        -----
        Since ``kwargs`` is a dictionary, the order of your arguments may not
        be preserved, and so the order of the new variables is not well
        defined. Assigning multiple variables within the same ``assign`` is
        possible, but you cannot reference other variables created within the
        same ``assign`` call.

        The new assigned variables that replace existing coordinates in the
        original dataset are still listed as coordinates in the returned
        Dataset.

        See Also
        --------
        pandas.DataFrame.assign

        Examples
        --------
        >>> x = xr.Dataset(
        ...     {
        ...         "temperature_c": (
        ...             ("lat", "lon"),
        ...             20 * np.random.rand(4).reshape(2, 2),
        ...         ),
        ...         "precipitation": (("lat", "lon"), np.random.rand(4).reshape(2, 2)),
        ...     },
        ...     coords={"lat": [10, 20], "lon": [150, 160]},
        ... )
        >>> x
        <xarray.Dataset> Size: 96B
        Dimensions:        (lat: 2, lon: 2)
        Coordinates:
          * lat            (lat) int64 16B 10 20
          * lon            (lon) int64 16B 150 160
        Data variables:
            temperature_c  (lat, lon) float64 32B 10.98 14.3 12.06 10.9
            precipitation  (lat, lon) float64 32B 0.4237 0.6459 0.4376 0.8918

        Where the value is a callable, evaluated on dataset:

        >>> x.assign(temperature_f=lambda x: x.temperature_c * 9 / 5 + 32)
        <xarray.Dataset> Size: 128B
        Dimensions:        (lat: 2, lon: 2)
        Coordinates:
          * lat            (lat) int64 16B 10 20
          * lon            (lon) int64 16B 150 160
        Data variables:
            temperature_c  (lat, lon) float64 32B 10.98 14.3 12.06 10.9
            precipitation  (lat, lon) float64 32B 0.4237 0.6459 0.4376 0.8918
            temperature_f  (lat, lon) float64 32B 51.76 57.75 53.7 51.62

        Alternatively, the same behavior can be achieved by directly referencing an existing dataarray:

        >>> x.assign(temperature_f=x["temperature_c"] * 9 / 5 + 32)
        <xarray.Dataset> Size: 128B
        Dimensions:        (lat: 2, lon: 2)
        Coordinates:
          * lat            (lat) int64 16B 10 20
          * lon            (lon) int64 16B 150 160
        Data variables:
            temperature_c  (lat, lon) float64 32B 10.98 14.3 12.06 10.9
            precipitation  (lat, lon) float64 32B 0.4237 0.6459 0.4376 0.8918
            temperature_f  (lat, lon) float64 32B 51.76 57.75 53.7 51.62

        """
        variables = either_dict_or_kwargs(variables, variables_kwargs, "assign")
        data = self.copy()

        # do all calculations first...
        results: CoercibleMapping = data._calc_assign_results(variables)

        # split data variables to add/replace vs. coordinates to replace
        results_data_vars: dict[Hashable, CoercibleValue] = {}
        results_coords: dict[Hashable, CoercibleValue] = {}
        for k, v in results.items():
            if k in data._coord_names:
                results_coords[k] = v
            else:
                results_data_vars[k] = v

        # ... and then assign
        data.coords.update(results_coords)
        data.update(results_data_vars)

        return data

    def to_dataarray(
        self, dim: Hashable = "variable", name: Hashable | None = None
    ) -> DataArray:
        """Convert this dataset into an xarray.DataArray

        The data variables of this dataset will be broadcast against each other
        and stacked along the first axis of the new array. All coordinates of
        this dataset will remain coordinates.

        Parameters
        ----------
        dim : Hashable, default: "variable"
            Name of the new dimension.
        name : Hashable or None, optional
            Name of the new data array.

        Returns
        -------
        array : xarray.DataArray
        """
        from xarray.core.dataarray import DataArray

        data_vars = [self.variables[k] for k in self.data_vars]
        broadcast_vars = broadcast_variables(*data_vars)
        data = duck_array_ops.stack([b.data for b in broadcast_vars], axis=0)

        dims = (dim,) + broadcast_vars[0].dims
        variable = Variable(dims, data, self.attrs, fastpath=True)

        coords = {k: v.variable for k, v in self.coords.items()}
        indexes = filter_indexes_from_coords(self._indexes, set(coords))
        new_dim_index = PandasIndex(list(self.data_vars), dim)
        indexes[dim] = new_dim_index
        coords.update(new_dim_index.create_variables())

        return DataArray._construct_direct(variable, coords, name, indexes)

    def to_array(
        self, dim: Hashable = "variable", name: Hashable | None = None
    ) -> DataArray:
        """Deprecated version of to_dataarray"""
        return self.to_dataarray(dim=dim, name=name)

    def _normalize_dim_order(
        self, dim_order: Sequence[Hashable] | None = None
    ) -> dict[Hashable, int]:
        """
        Check the validity of the provided dimensions if any and return the mapping
        between dimension name and their size.

        Parameters
        ----------
        dim_order: Sequence of Hashable or None, optional
            Dimension order to validate (default to the alphabetical order if None).

        Returns
        -------
        result : dict[Hashable, int]
            Validated dimensions mapping.

        """
        if dim_order is None:
            dim_order = list(self.dims)
        elif set(dim_order) != set(self.dims):
            raise ValueError(
                f"dim_order {dim_order} does not match the set of dimensions of this "
                f"Dataset: {list(self.dims)}"
            )

        ordered_dims = {k: self.sizes[k] for k in dim_order}

        return ordered_dims

    def to_pandas(self) -> pd.Series | pd.DataFrame:
        """Convert this dataset into a pandas object without changing the number of dimensions.

        The type of the returned object depends on the number of Dataset
        dimensions:

        * 0D -> `pandas.Series`
        * 1D -> `pandas.DataFrame`

        Only works for Datasets with 1 or fewer dimensions.
        """
        if len(self.dims) == 0:
            return pd.Series({k: v.item() for k, v in self.items()})
        if len(self.dims) == 1:
            return self.to_dataframe()
        raise ValueError(
            f"cannot convert Datasets with {len(self.dims)} dimensions into "
            "pandas objects without changing the number of dimensions. "
            "Please use Dataset.to_dataframe() instead."
        )

    def _to_dataframe(self, ordered_dims: Mapping[Any, int]):
        from xarray.core.extension_array import PandasExtensionArray

        columns_in_order = [k for k in self.variables if k not in self.dims]
        non_extension_array_columns = [
            k
            for k in columns_in_order
            if not pd.api.types.is_extension_array_dtype(self.variables[k].data)  # noqa: TID251
        ]
        extension_array_columns = [
            k
            for k in columns_in_order
            if pd.api.types.is_extension_array_dtype(self.variables[k].data)  # noqa: TID251
        ]
        extension_array_columns_different_index = [
            k
            for k in extension_array_columns
            if set(self.variables[k].dims) != set(ordered_dims.keys())
        ]
        extension_array_columns_same_index = [
            k
            for k in extension_array_columns
            if k not in extension_array_columns_different_index
        ]
        data = [
            self._variables[k].set_dims(ordered_dims).values.reshape(-1)
            for k in non_extension_array_columns
        ]
        index = self.coords.to_index([*ordered_dims])
        broadcasted_df = pd.DataFrame(
            {
                **dict(zip(non_extension_array_columns, data, strict=True)),
                **{
                    c: self.variables[c].data
                    for c in extension_array_columns_same_index
                },
            },
            index=index,
        )
        for extension_array_column in extension_array_columns_different_index:
            extension_array = self.variables[extension_array_column].data
            index = self[
                self.variables[extension_array_column].dims[0]
            ].coords.to_index()
            extension_array_df = pd.DataFrame(
                {extension_array_column: extension_array},
                index=pd.Index(index.array)
                if isinstance(index, PandasExtensionArray)  # type: ignore[redundant-expr]
                else index,
            )
            extension_array_df.index.name = self.variables[extension_array_column].dims[
                0
            ]
            broadcasted_df = broadcasted_df.join(extension_array_df)
        return broadcasted_df[columns_in_order]

    def to_dataframe(self, dim_order: Sequence[Hashable] | None = None) -> pd.DataFrame:
        """Convert this dataset into a pandas.DataFrame.

        Non-index variables in this dataset form the columns of the
        DataFrame. The DataFrame is indexed by the Cartesian product of
        this dataset's indices.

        Parameters
        ----------
        dim_order: Sequence of Hashable or None, optional
            Hierarchical dimension order for the resulting dataframe. All
            arrays are transposed to this order and then written out as flat
            vectors in contiguous order, so the last dimension in this list
            will be contiguous in the resulting DataFrame. This has a major
            influence on which operations are efficient on the resulting
            dataframe.

            If provided, must include all dimensions of this dataset. By
            default, dimensions are in the same order as in `Dataset.sizes`.

        Returns
        -------
        result : DataFrame
            Dataset as a pandas DataFrame.

        """

        ordered_dims = self._normalize_dim_order(dim_order=dim_order)

        return self._to_dataframe(ordered_dims=ordered_dims)

    def _set_sparse_data_from_dataframe(
        self, idx: pd.Index, arrays: list[tuple[Hashable, np.ndarray]], dims: tuple
    ) -> None:
        from sparse import COO

        if isinstance(idx, pd.MultiIndex):
            coords = np.stack([np.asarray(code) for code in idx.codes], axis=0)
            is_sorted = idx.is_monotonic_increasing
            shape = tuple(lev.size for lev in idx.levels)
        else:
            coords = np.arange(idx.size).reshape(1, -1)
            is_sorted = True
            shape = (idx.size,)

        for name, values in arrays:
            # In virtually all real use cases, the sparse array will now have
            # missing values and needs a fill_value. For consistency, don't
            # special case the rare exceptions (e.g., dtype=int without a
            # MultiIndex).
            dtype, fill_value = xrdtypes.maybe_promote(values.dtype)
            values = np.asarray(values, dtype=dtype)

            data = COO(
                coords,
                values,
                shape,
                has_duplicates=False,
                sorted=is_sorted,
                fill_value=fill_value,
            )
            self[name] = (dims, data)

    def _set_numpy_data_from_dataframe(
        self, idx: pd.Index, arrays: list[tuple[Hashable, np.ndarray]], dims: tuple
    ) -> None:
        if not isinstance(idx, pd.MultiIndex):
            for name, values in arrays:
                self[name] = (dims, values)
            return

        # NB: similar, more general logic, now exists in
        # variable.unstack_once; we could consider combining them at some
        # point.

        shape = tuple(lev.size for lev in idx.levels)
        indexer = tuple(idx.codes)

        # We already verified that the MultiIndex has all unique values, so
        # there are missing values if and only if the size of output arrays is
        # larger that the index.
        missing_values = math.prod(shape) > idx.shape[0]

        for name, values in arrays:
            # NumPy indexing is much faster than using DataFrame.reindex() to
            # fill in missing values:
            # https://stackoverflow.com/a/35049899/809705
            if missing_values:
                dtype, fill_value = xrdtypes.maybe_promote(values.dtype)
                data = np.full(shape, fill_value, dtype)
            else:
                # If there are no missing values, keep the existing dtype
                # instead of promoting to support NA, e.g., keep integer
                # columns as integers.
                # TODO: consider removing this special case, which doesn't
                # exist for sparse=True.
                data = np.zeros(shape, values.dtype)
            data[indexer] = values
            self[name] = (dims, data)

    @classmethod
    def from_dataframe(cls, dataframe: pd.DataFrame, sparse: bool = False) -> Self:
        """Convert a pandas.DataFrame into an xarray.Dataset

        Each column will be converted into an independent variable in the
        Dataset. If the dataframe's index is a MultiIndex, it will be expanded
        into a tensor product of one-dimensional indices (filling in missing
        values with NaN). If you rather preserve the MultiIndex use
        `xr.Dataset(df)`. This method will produce a Dataset very similar to
        that on which the 'to_dataframe' method was called, except with
        possibly redundant dimensions (since all dataset variables will have
        the same dimensionality).

        Parameters
        ----------
        dataframe : DataFrame
            DataFrame from which to copy data and indices.
        sparse : bool, default: False
            If true, create a sparse arrays instead of dense numpy arrays. This
            can potentially save a large amount of memory if the DataFrame has
            a MultiIndex. Requires the sparse package (sparse.pydata.org).

        Returns
        -------
        New Dataset.

        See Also
        --------
        xarray.DataArray.from_series
        pandas.DataFrame.to_xarray
        """
        # TODO: Add an option to remove dimensions along which the variables
        # are constant, to enable consistent serialization to/from a dataframe,
        # even if some variables have different dimensionality.

        if not dataframe.columns.is_unique:
            raise ValueError("cannot convert DataFrame with non-unique columns")

        idx = remove_unused_levels_categories(dataframe.index)

        if isinstance(idx, pd.MultiIndex) and not idx.is_unique:
            raise ValueError(
                "cannot convert a DataFrame with a non-unique MultiIndex into xarray"
            )

        arrays = []
        extension_arrays = []
        for k, v in dataframe.items():
            if not is_allowed_extension_array(v) or isinstance(
                v.array, UNSUPPORTED_EXTENSION_ARRAY_TYPES
            ):
                arrays.append((k, np.asarray(v)))
            else:
                extension_arrays.append((k, v))

        indexes: dict[Hashable, Index] = {}
        index_vars: dict[Hashable, Variable] = {}

        if isinstance(idx, pd.MultiIndex):
            dims = tuple(
                name if name is not None else f"level_{n}"  # type: ignore[redundant-expr,unused-ignore]
                for n, name in enumerate(idx.names)
            )
            for dim, lev in zip(dims, idx.levels, strict=True):
                xr_idx = PandasIndex(lev, dim)
                indexes[dim] = xr_idx
                index_vars.update(xr_idx.create_variables())
            arrays += [(k, np.asarray(v)) for k, v in extension_arrays]
            extension_arrays = []
        else:
            index_name = idx.name if idx.name is not None else "index"
            dims = (index_name,)
            xr_idx = PandasIndex(idx, index_name)
            indexes[index_name] = xr_idx
            index_vars.update(xr_idx.create_variables())

        obj = cls._construct_direct(index_vars, set(index_vars), indexes=indexes)

        if sparse:
            obj._set_sparse_data_from_dataframe(idx, arrays, dims)
        else:
            obj._set_numpy_data_from_dataframe(idx, arrays, dims)
        for name, extension_array in extension_arrays:
            obj[name] = (dims, extension_array)
        return obj[dataframe.columns] if len(dataframe.columns) else obj

    def to_dask_dataframe(
        self, dim_order: Sequence[Hashable] | None = None, set_index: bool = False
    ) -> DaskDataFrame:
        """
        Convert this dataset into a dask.dataframe.DataFrame.

        The dimensions, coordinates and data variables in this dataset form
        the columns of the DataFrame.

        Parameters
        ----------
        dim_order : list, optional
            Hierarchical dimension order for the resulting dataframe. All
            arrays are transposed to this order and then written out as flat
            vectors in contiguous order, so the last dimension in this list
            will be contiguous in the resulting DataFrame. This has a major
            influence on which operations are efficient on the resulting dask
            dataframe.

            If provided, must include all dimensions of this dataset. By
            default, dimensions are sorted alphabetically.
        set_index : bool, default: False
            If set_index=True, the dask DataFrame is indexed by this dataset's
            coordinate. Since dask DataFrames do not support multi-indexes,
            set_index only works if the dataset only contains one dimension.

        Returns
        -------
        dask.dataframe.DataFrame
        """

        import dask.array as da
        import dask.dataframe as dd

        ordered_dims = self._normalize_dim_order(dim_order=dim_order)

        columns = list(ordered_dims)
        columns.extend(k for k in self.coords if k not in self.dims)
        columns.extend(self.data_vars)

        ds_chunks = self.chunks

        series_list = []
        df_meta = pd.DataFrame()
        for name in columns:
            try:
                var = self.variables[name]
            except KeyError:
                # dimension without a matching coordinate
                size = self.sizes[name]
                data = da.arange(size, chunks=size, dtype=np.int64)
                var = Variable((name,), data)

            # IndexVariable objects have a dummy .chunk() method
            if isinstance(var, IndexVariable):
                var = var.to_base_variable()

            # Make sure var is a dask array, otherwise the array can become too large
            # when it is broadcasted to several dimensions:
            if not is_duck_dask_array(var._data):
                var = var.chunk()

            # Broadcast then flatten the array:
            var_new_dims = var.set_dims(ordered_dims).chunk(ds_chunks)
            dask_array = var_new_dims._data.reshape(-1)

            series = dd.from_dask_array(dask_array, columns=name, meta=df_meta)
            series_list.append(series)

        df = dd.concat(series_list, axis=1)

        if set_index:
            dim_order = [*ordered_dims]

            if len(dim_order) == 1:
                (dim,) = dim_order
                df = df.set_index(dim)
            else:
                # triggers an error about multi-indexes, even if only one
                # dimension is passed
                df = df.set_index(dim_order)

        return df

    def to_dict(
        self, data: bool | Literal["list", "array"] = "list", encoding: bool = False
    ) -> dict[str, Any]:
        """
        Convert this dataset to a dictionary following xarray naming
        conventions.

        Converts all variables and attributes to native Python objects
        Useful for converting to json. To avoid datetime incompatibility
        use decode_times=False kwarg in xarrray.open_dataset.

        Parameters
        ----------
        data : bool or {"list", "array"}, default: "list"
            Whether to include the actual data in the dictionary. When set to
            False, returns just the schema. If set to "array", returns data as
            underlying array type. If set to "list" (or True for backwards
            compatibility), returns data in lists of Python data types. Note
            that for obtaining the "list" output efficiently, use
            `ds.compute().to_dict(data="list")`.

        encoding : bool, default: False
            Whether to include the Dataset's encoding in the dictionary.

        Returns
        -------
        d : dict
            Dict with keys: "coords", "attrs", "dims", "data_vars" and optionally
            "encoding".

        See Also
        --------
        Dataset.from_dict
        DataArray.to_dict
        """
        d: dict = {
            "coords": {},
            "attrs": decode_numpy_dict_values(self.attrs),
            "dims": dict(self.sizes),
            "data_vars": {},
        }
        for k in self.coords:
            d["coords"].update(
                {k: self[k].variable.to_dict(data=data, encoding=encoding)}
            )
        for k in self.data_vars:
            d["data_vars"].update(
                {k: self[k].variable.to_dict(data=data, encoding=encoding)}
            )
        if encoding:
            d["encoding"] = dict(self.encoding)
        return d

    @classmethod
    def from_dict(cls, d: Mapping[Any, Any]) -> Self:
        """Convert a dictionary into an xarray.Dataset.

        Parameters
        ----------
        d : dict-like
            Mapping with a minimum structure of
                ``{"var_0": {"dims": [..], "data": [..]}, \
                            ...}``

        Returns
        -------
        obj : Dataset

        See also
        --------
        Dataset.to_dict
        DataArray.from_dict

        Examples
        --------
        >>> d = {
        ...     "t": {"dims": ("t"), "data": [0, 1, 2]},
        ...     "a": {"dims": ("t"), "data": ["a", "b", "c"]},
        ...     "b": {"dims": ("t"), "data": [10, 20, 30]},
        ... }
        >>> ds = xr.Dataset.from_dict(d)
        >>> ds
        <xarray.Dataset> Size: 60B
        Dimensions:  (t: 3)
        Coordinates:
          * t        (t) int64 24B 0 1 2
        Data variables:
            a        (t) <U1 12B 'a' 'b' 'c'
            b        (t) int64 24B 10 20 30

        >>> d = {
        ...     "coords": {
        ...         "t": {"dims": "t", "data": [0, 1, 2], "attrs": {"units": "s"}}
        ...     },
        ...     "attrs": {"title": "air temperature"},
        ...     "dims": "t",
        ...     "data_vars": {
        ...         "a": {"dims": "t", "data": [10, 20, 30]},
        ...         "b": {"dims": "t", "data": ["a", "b", "c"]},
        ...     },
        ... }
        >>> ds = xr.Dataset.from_dict(d)
        >>> ds
        <xarray.Dataset> Size: 60B
        Dimensions:  (t: 3)
        Coordinates:
          * t        (t) int64 24B 0 1 2
        Data variables:
            a        (t) int64 24B 10 20 30
            b        (t) <U1 12B 'a' 'b' 'c'
        Attributes:
            title:    air temperature

        """

        variables: Iterable[tuple[Hashable, Any]]
        if not {"coords", "data_vars"}.issubset(set(d)):
            variables = d.items()
        else:
            import itertools

            variables = itertools.chain(
                d.get("coords", {}).items(), d.get("data_vars", {}).items()
            )
        try:
            variable_dict = {
                k: (v["dims"], v["data"], v.get("attrs"), v.get("encoding"))
                for k, v in variables
            }
        except KeyError as e:
            raise ValueError(
                f"cannot convert dict without the key '{e.args[0]}'"
            ) from e
        obj = cls(variable_dict)

        # what if coords aren't dims?
        coords = set(d.get("coords", {})) - set(d.get("dims", {}))
        obj = obj.set_coords(coords)

        obj.attrs.update(d.get("attrs", {}))
        obj.encoding.update(d.get("encoding", {}))

        return obj

    def _unary_op(self, f, *args, **kwargs) -> Self:
        variables = {}
        keep_attrs = kwargs.pop("keep_attrs", None)
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)
        for k, v in self._variables.items():
            if k in self._coord_names:
                variables[k] = v
            else:
                variables[k] = f(v, *args, **kwargs)
                if keep_attrs:
                    variables[k]._attrs = v._attrs
        attrs = self._attrs if keep_attrs else None
        return self._replace_with_new_dims(variables, attrs=attrs)

    def _binary_op(self, other, f, reflexive=False, join=None) -> Dataset:
        from xarray.core.dataarray import DataArray
        from xarray.core.datatree import DataTree
        from xarray.core.groupby import GroupBy

        if isinstance(other, DataTree | GroupBy):
            return NotImplemented
        align_type = OPTIONS["arithmetic_join"] if join is None else join
        if isinstance(other, DataArray | Dataset):
            self, other = align(self, other, join=align_type, copy=False)
        g = f if not reflexive else lambda x, y: f(y, x)
        ds = self._calculate_binary_op(g, other, join=align_type)
        keep_attrs = _get_keep_attrs(default=False)
        if keep_attrs:
            ds.attrs = self.attrs
        return ds

    def _inplace_binary_op(self, other, f) -> Self:
        from xarray.core.dataarray import DataArray
        from xarray.core.groupby import GroupBy

        if isinstance(other, GroupBy):
            raise TypeError(
                "in-place operations between a Dataset and "
                "a grouped object are not permitted"
            )
        # we don't actually modify arrays in-place with in-place Dataset
        # arithmetic -- this lets us automatically align things
        if isinstance(other, DataArray | Dataset):
            other = other.reindex_like(self, copy=False)
        g = ops.inplace_to_noninplace_op(f)
        ds = self._calculate_binary_op(g, other, inplace=True)
        self._replace_with_new_dims(
            ds._variables,
            ds._coord_names,
            attrs=ds._attrs,
            indexes=ds._indexes,
            inplace=True,
        )
        return self

    def _calculate_binary_op(
        self, f, other, join="inner", inplace: bool = False
    ) -> Dataset:
        def apply_over_both(lhs_data_vars, rhs_data_vars, lhs_vars, rhs_vars):
            if inplace and set(lhs_data_vars) != set(rhs_data_vars):
                raise ValueError(
                    "datasets must have the same data variables "
                    f"for in-place arithmetic operations: {list(lhs_data_vars)}, {list(rhs_data_vars)}"
                )

            dest_vars = {}

            for k in lhs_data_vars:
                if k in rhs_data_vars:
                    dest_vars[k] = f(lhs_vars[k], rhs_vars[k])
                elif join in ["left", "outer"]:
                    dest_vars[k] = f(lhs_vars[k], np.nan)
            for k in rhs_data_vars:
                if k not in dest_vars and join in ["right", "outer"]:
                    dest_vars[k] = f(rhs_vars[k], np.nan)
            return dest_vars

        if utils.is_dict_like(other) and not isinstance(other, Dataset):
            # can't use our shortcut of doing the binary operation with
            # Variable objects, so apply over our data vars instead.
            new_data_vars = apply_over_both(
                self.data_vars, other, self.data_vars, other
            )
            return type(self)(new_data_vars)

        other_coords: Coordinates | None = getattr(other, "coords", None)
        ds = self.coords.merge(other_coords)

        if isinstance(other, Dataset):
            new_vars = apply_over_both(
                self.data_vars, other.data_vars, self.variables, other.variables
            )
        else:
            other_variable = getattr(other, "variable", other)
            new_vars = {k: f(self.variables[k], other_variable) for k in self.data_vars}
        ds._variables.update(new_vars)
        ds._dims = calculate_dimensions(ds._variables)
        return ds

    def _copy_attrs_from(self, other):
        self.attrs = other.attrs
        for v in other.variables:
            if v in self.variables:
                self.variables[v].attrs = other.variables[v].attrs

    def diff(
        self,
        dim: Hashable,
        n: int = 1,
        *,
        label: Literal["upper", "lower"] = "upper",
    ) -> Self:
        """Calculate the n-th order discrete difference along given axis.

        Parameters
        ----------
        dim : Hashable
            Dimension over which to calculate the finite difference.
        n : int, default: 1
            The number of times values are differenced.
        label : {"upper", "lower"}, default: "upper"
            The new coordinate in dimension ``dim`` will have the
            values of either the minuend's or subtrahend's coordinate
            for values 'upper' and 'lower', respectively.

        Returns
        -------
        difference : Dataset
            The n-th order finite difference of this object.

        Notes
        -----
        `n` matches numpy's behavior and is different from pandas' first argument named
        `periods`.

        Examples
        --------
        >>> ds = xr.Dataset({"foo": ("x", [5, 5, 6, 6])})
        >>> ds.diff("x")
        <xarray.Dataset> Size: 24B
        Dimensions:  (x: 3)
        Dimensions without coordinates: x
        Data variables:
            foo      (x) int64 24B 0 1 0
        >>> ds.diff("x", 2)
        <xarray.Dataset> Size: 16B
        Dimensions:  (x: 2)
        Dimensions without coordinates: x
        Data variables:
            foo      (x) int64 16B 1 -1

        See Also
        --------
        Dataset.differentiate
        """
        if n == 0:
            return self
        if n < 0:
            raise ValueError(f"order `n` must be non-negative but got {n}")

        # prepare slices
        slice_start = {dim: slice(None, -1)}
        slice_end = {dim: slice(1, None)}

        # prepare new coordinate
        if label == "upper":
            slice_new = slice_end
        elif label == "lower":
            slice_new = slice_start
        else:
            raise ValueError("The 'label' argument has to be either 'upper' or 'lower'")

        indexes, index_vars = isel_indexes(self.xindexes, slice_new)
        variables = {}

        for name, var in self.variables.items():
            if name in index_vars:
                variables[name] = index_vars[name]
            elif dim in var.dims:
                if name in self.data_vars:
                    variables[name] = var.isel(slice_end) - var.isel(slice_start)
                else:
                    variables[name] = var.isel(slice_new)
            else:
                variables[name] = var

        difference = self._replace_with_new_dims(variables, indexes=indexes)

        if n > 1:
            return difference.diff(dim, n - 1)
        else:
            return difference

    def shift(
        self,
        shifts: Mapping[Any, int] | None = None,
        fill_value: Any = xrdtypes.NA,
        **shifts_kwargs: int,
    ) -> Self:
        """Shift this dataset by an offset along one or more dimensions.

        Only data variables are moved; coordinates stay in place. This is
        consistent with the behavior of ``shift`` in pandas.

        Values shifted from beyond array bounds will appear at one end of
        each dimension, which are filled according to `fill_value`. For periodic
        offsets instead see `roll`.

        Parameters
        ----------
        shifts : mapping of hashable to int
            Integer offset to shift along each of the given dimensions.
            Positive offsets shift to the right; negative offsets shift to the
            left.
        fill_value : scalar or dict-like, optional
            Value to use for newly missing values. If a dict-like, maps
            variable names (including coordinates) to fill values.
        **shifts_kwargs
            The keyword arguments form of ``shifts``.
            One of shifts or shifts_kwargs must be provided.

        Returns
        -------
        shifted : Dataset
            Dataset with the same coordinates and attributes but shifted data
            variables.

        See Also
        --------
        roll

        Examples
        --------
        >>> ds = xr.Dataset({"foo": ("x", list("abcde"))})
        >>> ds.shift(x=2)
        <xarray.Dataset> Size: 40B
        Dimensions:  (x: 5)
        Dimensions without coordinates: x
        Data variables:
            foo      (x) object 40B nan nan 'a' 'b' 'c'
        """
        shifts = either_dict_or_kwargs(shifts, shifts_kwargs, "shift")
        invalid = tuple(k for k in shifts if k not in self.dims)
        if invalid:
            raise ValueError(
                f"Dimensions {invalid} not found in data dimensions {tuple(self.dims)}"
            )

        variables = {}
        for name, var in self.variables.items():
            if name in self.data_vars:
                fill_value_ = (
                    fill_value.get(name, xrdtypes.NA)
                    if isinstance(fill_value, dict)
                    else fill_value
                )

                var_shifts = {k: v for k, v in shifts.items() if k in var.dims}
                variables[name] = var.shift(fill_value=fill_value_, shifts=var_shifts)
            else:
                variables[name] = var

        return self._replace(variables)

    def roll(
        self,
        shifts: Mapping[Any, int] | None = None,
        roll_coords: bool = False,
        **shifts_kwargs: int,
    ) -> Self:
        """Roll this dataset by an offset along one or more dimensions.

        Unlike shift, roll treats the given dimensions as periodic, so will not
        create any missing values to be filled.

        Also unlike shift, roll may rotate all variables, including coordinates
        if specified. The direction of rotation is consistent with
        :py:func:`numpy.roll`.

        Parameters
        ----------
        shifts : mapping of hashable to int, optional
            A dict with keys matching dimensions and values given
            by integers to rotate each of the given dimensions. Positive
            offsets roll to the right; negative offsets roll to the left.
        roll_coords : bool, default: False
            Indicates whether to roll the coordinates by the offset too.
        **shifts_kwargs : {dim: offset, ...}, optional
            The keyword arguments form of ``shifts``.
            One of shifts or shifts_kwargs must be provided.

        Returns
        -------
        rolled : Dataset
            Dataset with the same attributes but rolled data and coordinates.

        See Also
        --------
        shift

        Examples
        --------
        >>> ds = xr.Dataset({"foo": ("x", list("abcde"))}, coords={"x": np.arange(5)})
        >>> ds.roll(x=2)
        <xarray.Dataset> Size: 60B
        Dimensions:  (x: 5)
        Coordinates:
          * x        (x) int64 40B 0 1 2 3 4
        Data variables:
            foo      (x) <U1 20B 'd' 'e' 'a' 'b' 'c'

        >>> ds.roll(x=2, roll_coords=True)
        <xarray.Dataset> Size: 60B
        Dimensions:  (x: 5)
        Coordinates:
          * x        (x) int64 40B 3 4 0 1 2
        Data variables:
            foo      (x) <U1 20B 'd' 'e' 'a' 'b' 'c'

        """
        shifts = either_dict_or_kwargs(shifts, shifts_kwargs, "roll")
        invalid = [k for k in shifts if k not in self.dims]
        if invalid:
            raise ValueError(
                f"Dimensions {invalid} not found in data dimensions {tuple(self.dims)}"
            )

        unrolled_vars: tuple[Hashable, ...]

        if roll_coords:
            indexes, index_vars = roll_indexes(self.xindexes, shifts)
            unrolled_vars = ()
        else:
            indexes = dict(self._indexes)
            index_vars = dict(self.xindexes.variables)
            unrolled_vars = tuple(self.coords)

        variables = {}
        for k, var in self.variables.items():
            if k in index_vars:
                variables[k] = index_vars[k]
            elif k not in unrolled_vars:
                variables[k] = var.roll(
                    shifts={k: s for k, s in shifts.items() if k in var.dims}
                )
            else:
                variables[k] = var

        return self._replace(variables, indexes=indexes)

    def sortby(
        self,
        variables: (
            Hashable
            | DataArray
            | Sequence[Hashable | DataArray]
            | Callable[[Self], Hashable | DataArray | list[Hashable | DataArray]]
        ),
        ascending: bool = True,
    ) -> Self:
        """
        Sort object by labels or values (along an axis).

        Sorts the dataset, either along specified dimensions,
        or according to values of 1-D dataarrays that share dimension
        with calling object.

        If the input variables are dataarrays, then the dataarrays are aligned
        (via left-join) to the calling object prior to sorting by cell values.
        NaNs are sorted to the end, following Numpy convention.

        If multiple sorts along the same dimension is
        given, numpy's lexsort is performed along that dimension:
        https://numpy.org/doc/stable/reference/generated/numpy.lexsort.html
        and the FIRST key in the sequence is used as the primary sort key,
        followed by the 2nd key, etc.

        Parameters
        ----------
        variables : Hashable, DataArray, sequence of Hashable or DataArray, or Callable
            1D DataArray objects or name(s) of 1D variable(s) in coords whose values are
            used to sort this array. If a callable, the callable is passed this object,
            and the result is used as the value for cond.
        ascending : bool, default: True
            Whether to sort by ascending or descending order.

        Returns
        -------
        sorted : Dataset
            A new dataset where all the specified dims are sorted by dim
            labels.

        See Also
        --------
        DataArray.sortby
        numpy.sort
        pandas.sort_values
        pandas.sort_index

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     {
        ...         "A": (("x", "y"), [[1, 2], [3, 4]]),
        ...         "B": (("x", "y"), [[5, 6], [7, 8]]),
        ...     },
        ...     coords={"x": ["b", "a"], "y": [1, 0]},
        ... )
        >>> ds.sortby("x")
        <xarray.Dataset> Size: 88B
        Dimensions:  (x: 2, y: 2)
        Coordinates:
          * x        (x) <U1 8B 'a' 'b'
          * y        (y) int64 16B 1 0
        Data variables:
            A        (x, y) int64 32B 3 4 1 2
            B        (x, y) int64 32B 7 8 5 6
        >>> ds.sortby(lambda x: -x["y"])
        <xarray.Dataset> Size: 88B
        Dimensions:  (x: 2, y: 2)
        Coordinates:
          * x        (x) <U1 8B 'b' 'a'
          * y        (y) int64 16B 1 0
        Data variables:
            A        (x, y) int64 32B 1 2 3 4
            B        (x, y) int64 32B 5 6 7 8
        """
        from xarray.core.dataarray import DataArray

        if callable(variables):
            variables = variables(self)
        if not isinstance(variables, list):
            variables = [variables]
        arrays = [v if isinstance(v, DataArray) else self[v] for v in variables]
        aligned_vars = align(self, *arrays, join="left")
        aligned_self = cast("Self", aligned_vars[0])
        aligned_other_vars = cast(tuple[DataArray, ...], aligned_vars[1:])
        vars_by_dim = defaultdict(list)
        for data_array in aligned_other_vars:
            if data_array.ndim != 1:
                raise ValueError("Input DataArray is not 1-D.")
            (key,) = data_array.dims
            vars_by_dim[key].append(data_array)

        indices = {}
        for key, arrays in vars_by_dim.items():
            order = np.lexsort(tuple(reversed(arrays)))
            indices[key] = order if ascending else order[::-1]
        return aligned_self.isel(indices)

    def quantile(
        self,
        q: ArrayLike,
        dim: Dims = None,
        *,
        method: QuantileMethods = "linear",
        numeric_only: bool = False,
        keep_attrs: bool | None = None,
        skipna: bool | None = None,
        interpolation: QuantileMethods | None = None,
    ) -> Self:
        """Compute the qth quantile of the data along the specified dimension.

        Returns the qth quantiles(s) of the array elements for each variable
        in the Dataset.

        Parameters
        ----------
        q : float or array-like of float
            Quantile to compute, which must be between 0 and 1 inclusive.
        dim : str or Iterable of Hashable, optional
            Dimension(s) over which to apply quantile.
        method : str, default: "linear"
            This optional parameter specifies the interpolation method to use when the
            desired quantile lies between two data points. The options sorted by their R
            type as summarized in the H&F paper [1]_ are:

                1. "inverted_cdf"
                2. "averaged_inverted_cdf"
                3. "closest_observation"
                4. "interpolated_inverted_cdf"
                5. "hazen"
                6. "weibull"
                7. "linear"  (default)
                8. "median_unbiased"
                9. "normal_unbiased"

            The first three methods are discontiuous.  The following discontinuous
            variations of the default "linear" (7.) option are also available:

                * "lower"
                * "higher"
                * "midpoint"
                * "nearest"

            See :py:func:`numpy.quantile` or [1]_ for details. The "method" argument
            was previously called "interpolation", renamed in accordance with numpy
            version 1.22.0.

        keep_attrs : bool, optional
            If True, the dataset's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        numeric_only : bool, optional
            If True, only apply ``func`` to variables with a numeric dtype.
        skipna : bool, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or skipna=True has not been
            implemented (object, datetime64 or timedelta64).

        Returns
        -------
        quantiles : Dataset
            If `q` is a single quantile, then the result is a scalar for each
            variable in data_vars. If multiple percentiles are given, first
            axis of the result corresponds to the quantile and a quantile
            dimension is added to the return Dataset. The other dimensions are
            the dimensions that remain after the reduction of the array.

        See Also
        --------
        numpy.nanquantile, numpy.quantile, pandas.Series.quantile, DataArray.quantile

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     {"a": (("x", "y"), [[0.7, 4.2, 9.4, 1.5], [6.5, 7.3, 2.6, 1.9]])},
        ...     coords={"x": [7, 9], "y": [1, 1.5, 2, 2.5]},
        ... )
        >>> ds.quantile(0)  # or ds.quantile(0, dim=...)
        <xarray.Dataset> Size: 16B
        Dimensions:   ()
        Coordinates:
            quantile  float64 8B 0.0
        Data variables:
            a         float64 8B 0.7
        >>> ds.quantile(0, dim="x")
        <xarray.Dataset> Size: 72B
        Dimensions:   (y: 4)
        Coordinates:
          * y         (y) float64 32B 1.0 1.5 2.0 2.5
            quantile  float64 8B 0.0
        Data variables:
            a         (y) float64 32B 0.7 4.2 2.6 1.5
        >>> ds.quantile([0, 0.5, 1])
        <xarray.Dataset> Size: 48B
        Dimensions:   (quantile: 3)
        Coordinates:
          * quantile  (quantile) float64 24B 0.0 0.5 1.0
        Data variables:
            a         (quantile) float64 24B 0.7 3.4 9.4
        >>> ds.quantile([0, 0.5, 1], dim="x")
        <xarray.Dataset> Size: 152B
        Dimensions:   (quantile: 3, y: 4)
        Coordinates:
          * y         (y) float64 32B 1.0 1.5 2.0 2.5
          * quantile  (quantile) float64 24B 0.0 0.5 1.0
        Data variables:
            a         (quantile, y) float64 96B 0.7 4.2 2.6 1.5 3.6 ... 6.5 7.3 9.4 1.9

        References
        ----------
        .. [1] R. J. Hyndman and Y. Fan,
           "Sample quantiles in statistical packages,"
           The American Statistician, 50(4), pp. 361-365, 1996
        """

        # interpolation renamed to method in version 0.21.0
        # check here and in variable to avoid repeated warnings
        if interpolation is not None:
            warnings.warn(
                "The `interpolation` argument to quantile was renamed to `method`.",
                FutureWarning,
                stacklevel=2,
            )

            if method != "linear":
                raise TypeError("Cannot pass interpolation and method keywords!")

            method = interpolation

        dims: set[Hashable]
        if isinstance(dim, str):
            dims = {dim}
        elif dim is None or dim is ...:
            dims = set(self.dims)
        else:
            dims = set(dim)

        invalid_dims = set(dims) - set(self.dims)
        if invalid_dims:
            raise ValueError(
                f"Dimensions {tuple(invalid_dims)} not found in data dimensions {tuple(self.dims)}"
            )

        q = np.asarray(q, dtype=np.float64)

        variables = {}
        for name, var in self.variables.items():
            reduce_dims = [d for d in var.dims if d in dims]
            if reduce_dims or not var.dims:
                if name not in self.coords and (
                    not numeric_only
                    or np.issubdtype(var.dtype, np.number)
                    or var.dtype == np.bool_
                ):
                    variables[name] = var.quantile(
                        q,
                        dim=reduce_dims,
                        method=method,
                        keep_attrs=keep_attrs,
                        skipna=skipna,
                    )

            else:
                variables[name] = var

        # construct the new dataset
        coord_names = {k for k in self.coords if k in variables}
        indexes = {k: v for k, v in self._indexes.items() if k in variables}
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)
        attrs = self.attrs if keep_attrs else None
        new = self._replace_with_new_dims(
            variables, coord_names=coord_names, attrs=attrs, indexes=indexes
        )
        return new.assign_coords(quantile=q)

    def rank(
        self,
        dim: Hashable,
        *,
        pct: bool = False,
        keep_attrs: bool | None = None,
    ) -> Self:
        """Ranks the data.

        Equal values are assigned a rank that is the average of the ranks that
        would have been otherwise assigned to all of the values within
        that set.
        Ranks begin at 1, not 0. If pct is True, computes percentage ranks.

        NaNs in the input array are returned as NaNs.

        The `bottleneck` library is required.

        Parameters
        ----------
        dim : Hashable
            Dimension over which to compute rank.
        pct : bool, default: False
            If True, compute percentage ranks, otherwise compute integer ranks.
        keep_attrs : bool or None, optional
            If True, the dataset's attributes (`attrs`) will be copied from
            the original object to the new one.  If False, the new
            object will be returned without attributes.

        Returns
        -------
        ranked : Dataset
            Variables that do not depend on `dim` are dropped.
        """
        if not OPTIONS["use_bottleneck"]:
            raise RuntimeError(
                "rank requires bottleneck to be enabled."
                " Call `xr.set_options(use_bottleneck=True)` to enable it."
            )

        if dim not in self.dims:
            raise ValueError(
                f"Dimension {dim!r} not found in data dimensions {tuple(self.dims)}"
            )

        variables = {}
        for name, var in self.variables.items():
            if name in self.data_vars:
                if dim in var.dims:
                    variables[name] = var.rank(dim, pct=pct)
            else:
                variables[name] = var

        coord_names = set(self.coords)
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)
        attrs = self.attrs if keep_attrs else None
        return self._replace(variables, coord_names, attrs=attrs)

    def differentiate(
        self,
        coord: Hashable,
        edge_order: Literal[1, 2] = 1,
        datetime_unit: DatetimeUnitOptions | None = None,
    ) -> Self:
        """Differentiate with the second order accurate central
        differences.

        .. note::
            This feature is limited to simple cartesian geometry, i.e. coord
            must be one dimensional.

        Parameters
        ----------
        coord : Hashable
            The coordinate to be used to compute the gradient.
        edge_order : {1, 2}, default: 1
            N-th order accurate differences at the boundaries.
        datetime_unit : None or {"W", "D", "h", "m", "s", "ms", \
            "us", "ns", "ps", "fs", "as", None}, default: None
            Unit to compute gradient. Only valid for datetime coordinate.

        Returns
        -------
        differentiated: Dataset

        See also
        --------
        numpy.gradient: corresponding numpy function
        """
        if coord not in self.variables and coord not in self.dims:
            variables_and_dims = tuple(set(self.variables.keys()).union(self.dims))
            raise ValueError(
                f"Coordinate {coord!r} not found in variables or dimensions {variables_and_dims}."
            )

        coord_var = self[coord].variable
        if coord_var.ndim != 1:
            raise ValueError(
                f"Coordinate {coord} must be 1 dimensional but is {coord_var.ndim}"
                " dimensional"
            )

        dim = coord_var.dims[0]
        if _contains_datetime_like_objects(coord_var):
            if coord_var.dtype.kind in "mM" and datetime_unit is None:
                datetime_unit = cast(
                    "DatetimeUnitOptions", np.datetime_data(coord_var.dtype)[0]
                )
            elif datetime_unit is None:
                datetime_unit = "s"  # Default to seconds for cftime objects
            coord_var = coord_var._to_numeric(datetime_unit=datetime_unit)

        variables = {}
        for k, v in self.variables.items():
            if k in self.data_vars and dim in v.dims and k not in self.coords:
                if _contains_datetime_like_objects(v):
                    v = v._to_numeric(datetime_unit=datetime_unit)
                grad = duck_array_ops.gradient(
                    v.data,
                    coord_var.data,
                    edge_order=edge_order,
                    axis=v.get_axis_num(dim),
                )
                variables[k] = Variable(v.dims, grad)
            else:
                variables[k] = v
        return self._replace(variables)

    def integrate(
        self,
        coord: Hashable | Sequence[Hashable],
        datetime_unit: DatetimeUnitOptions = None,
    ) -> Self:
        """Integrate along the given coordinate using the trapezoidal rule.

        .. note::
            This feature is limited to simple cartesian geometry, i.e. coord
            must be one dimensional.

        Parameters
        ----------
        coord : hashable, or sequence of hashable
            Coordinate(s) used for the integration.
        datetime_unit : {'W', 'D', 'h', 'm', 's', 'ms', 'us', 'ns', \
                        'ps', 'fs', 'as', None}, optional
            Specify the unit if datetime coordinate is used.

        Returns
        -------
        integrated : Dataset

        See also
        --------
        DataArray.integrate
        numpy.trapz : corresponding numpy function

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     data_vars={"a": ("x", [5, 5, 6, 6]), "b": ("x", [1, 2, 1, 0])},
        ...     coords={"x": [0, 1, 2, 3], "y": ("x", [1, 7, 3, 5])},
        ... )
        >>> ds
        <xarray.Dataset> Size: 128B
        Dimensions:  (x: 4)
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
            y        (x) int64 32B 1 7 3 5
        Data variables:
            a        (x) int64 32B 5 5 6 6
            b        (x) int64 32B 1 2 1 0
        >>> ds.integrate("x")
        <xarray.Dataset> Size: 16B
        Dimensions:  ()
        Data variables:
            a        float64 8B 16.5
            b        float64 8B 3.5
        >>> ds.integrate("y")
        <xarray.Dataset> Size: 16B
        Dimensions:  ()
        Data variables:
            a        float64 8B 20.0
            b        float64 8B 4.0
        """
        if not isinstance(coord, list | tuple):
            coord = (coord,)
        result = self
        for c in coord:
            result = result._integrate_one(c, datetime_unit=datetime_unit)
        return result

    def _integrate_one(self, coord, datetime_unit=None, cumulative=False):
        from xarray.core.variable import Variable

        if coord not in self.variables and coord not in self.dims:
            variables_and_dims = tuple(set(self.variables.keys()).union(self.dims))
            raise ValueError(
                f"Coordinate {coord!r} not found in variables or dimensions {variables_and_dims}."
            )

        coord_var = self[coord].variable
        if coord_var.ndim != 1:
            raise ValueError(
                f"Coordinate {coord} must be 1 dimensional but is {coord_var.ndim}"
                " dimensional"
            )

        dim = coord_var.dims[0]
        if _contains_datetime_like_objects(coord_var):
            if coord_var.dtype.kind in "mM" and datetime_unit is None:
                datetime_unit, _ = np.datetime_data(coord_var.dtype)
            elif datetime_unit is None:
                datetime_unit = "s"  # Default to seconds for cftime objects
            coord_var = coord_var._replace(
                data=datetime_to_numeric(coord_var.data, datetime_unit=datetime_unit)
            )

        variables = {}
        coord_names = set()
        for k, v in self.variables.items():
            if k in self.coords:
                if dim not in v.dims or cumulative:
                    variables[k] = v
                    coord_names.add(k)
            elif k in self.data_vars and dim in v.dims:
                coord_data = to_like_array(coord_var.data, like=v.data)
                if _contains_datetime_like_objects(v):
                    v = datetime_to_numeric(v, datetime_unit=datetime_unit)
                if cumulative:
                    integ = duck_array_ops.cumulative_trapezoid(
                        v.data, coord_data, axis=v.get_axis_num(dim)
                    )
                    v_dims = v.dims
                else:
                    integ = duck_array_ops.trapz(
                        v.data, coord_data, axis=v.get_axis_num(dim)
                    )
                    v_dims = list(v.dims)
                    v_dims.remove(dim)
                variables[k] = Variable(v_dims, integ)
            else:
                variables[k] = v
        indexes = {k: v for k, v in self._indexes.items() if k in variables}
        return self._replace_with_new_dims(
            variables, coord_names=coord_names, indexes=indexes
        )

    def cumulative_integrate(
        self,
        coord: Hashable | Sequence[Hashable],
        datetime_unit: DatetimeUnitOptions = None,
    ) -> Self:
        """Integrate along the given coordinate using the trapezoidal rule.

        .. note::
            This feature is limited to simple cartesian geometry, i.e. coord
            must be one dimensional.

            The first entry of the cumulative integral of each variable is always 0, in
            order to keep the length of the dimension unchanged between input and
            output.

        Parameters
        ----------
        coord : hashable, or sequence of hashable
            Coordinate(s) used for the integration.
        datetime_unit : {'W', 'D', 'h', 'm', 's', 'ms', 'us', 'ns', \
                        'ps', 'fs', 'as', None}, optional
            Specify the unit if datetime coordinate is used.

        Returns
        -------
        integrated : Dataset

        See also
        --------
        DataArray.cumulative_integrate
        scipy.integrate.cumulative_trapezoid : corresponding scipy function

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     data_vars={"a": ("x", [5, 5, 6, 6]), "b": ("x", [1, 2, 1, 0])},
        ...     coords={"x": [0, 1, 2, 3], "y": ("x", [1, 7, 3, 5])},
        ... )
        >>> ds
        <xarray.Dataset> Size: 128B
        Dimensions:  (x: 4)
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
            y        (x) int64 32B 1 7 3 5
        Data variables:
            a        (x) int64 32B 5 5 6 6
            b        (x) int64 32B 1 2 1 0
        >>> ds.cumulative_integrate("x")
        <xarray.Dataset> Size: 128B
        Dimensions:  (x: 4)
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
            y        (x) int64 32B 1 7 3 5
        Data variables:
            a        (x) float64 32B 0.0 5.0 10.5 16.5
            b        (x) float64 32B 0.0 1.5 3.0 3.5
        >>> ds.cumulative_integrate("y")
        <xarray.Dataset> Size: 128B
        Dimensions:  (x: 4)
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
            y        (x) int64 32B 1 7 3 5
        Data variables:
            a        (x) float64 32B 0.0 30.0 8.0 20.0
            b        (x) float64 32B 0.0 9.0 3.0 4.0
        """
        if not isinstance(coord, list | tuple):
            coord = (coord,)
        result = self
        for c in coord:
            result = result._integrate_one(
                c, datetime_unit=datetime_unit, cumulative=True
            )
        return result

    @property
    def real(self) -> Self:
        """
        The real part of each data variable.

        See Also
        --------
        numpy.ndarray.real
        """
        return self.map(lambda x: x.real, keep_attrs=True)

    @property
    def imag(self) -> Self:
        """
        The imaginary part of each data variable.

        See Also
        --------
        numpy.ndarray.imag
        """
        return self.map(lambda x: x.imag, keep_attrs=True)

    plot = utils.UncachedAccessor(DatasetPlotAccessor)

    def filter_by_attrs(self, **kwargs) -> Self:
        """Returns a ``Dataset`` with variables that match specific conditions.

        Can pass in ``key=value`` or ``key=callable``.  A Dataset is returned
        containing only the variables for which all the filter tests pass.
        These tests are either ``key=value`` for which the attribute ``key``
        has the exact value ``value`` or the callable passed into
        ``key=callable`` returns True. The callable will be passed a single
        value, either the value of the attribute ``key`` or ``None`` if the
        DataArray does not have an attribute with the name ``key``.

        Parameters
        ----------
        **kwargs
            key : str
                Attribute name.
            value : callable or obj
                If value is a callable, it should return a boolean in the form
                of bool = func(attr) where attr is da.attrs[key].
                Otherwise, value will be compared to the each
                DataArray's attrs[key].

        Returns
        -------
        new : Dataset
            New dataset with variables filtered by attribute.

        Examples
        --------
        >>> temp = 15 + 8 * np.random.randn(2, 2, 3)
        >>> precip = 10 * np.random.rand(2, 2, 3)
        >>> lon = [[-99.83, -99.32], [-99.79, -99.23]]
        >>> lat = [[42.25, 42.21], [42.63, 42.59]]
        >>> dims = ["x", "y", "time"]
        >>> temp_attr = dict(standard_name="air_potential_temperature")
        >>> precip_attr = dict(standard_name="convective_precipitation_flux")

        >>> ds = xr.Dataset(
        ...     dict(
        ...         temperature=(dims, temp, temp_attr),
        ...         precipitation=(dims, precip, precip_attr),
        ...     ),
        ...     coords=dict(
        ...         lon=(["x", "y"], lon),
        ...         lat=(["x", "y"], lat),
        ...         time=pd.date_range("2014-09-06", periods=3),
        ...         reference_time=pd.Timestamp("2014-09-05"),
        ...     ),
        ... )

        Get variables matching a specific standard_name:

        >>> ds.filter_by_attrs(standard_name="convective_precipitation_flux")
        <xarray.Dataset> Size: 192B
        Dimensions:         (x: 2, y: 2, time: 3)
        Coordinates:
            lon             (x, y) float64 32B -99.83 -99.32 -99.79 -99.23
            lat             (x, y) float64 32B 42.25 42.21 42.63 42.59
          * time            (time) datetime64[ns] 24B 2014-09-06 2014-09-07 2014-09-08
            reference_time  datetime64[ns] 8B 2014-09-05
        Dimensions without coordinates: x, y
        Data variables:
            precipitation   (x, y, time) float64 96B 5.68 9.256 0.7104 ... 4.615 7.805

        Get all variables that have a standard_name attribute:

        >>> standard_name = lambda v: v is not None
        >>> ds.filter_by_attrs(standard_name=standard_name)
        <xarray.Dataset> Size: 288B
        Dimensions:         (x: 2, y: 2, time: 3)
        Coordinates:
            lon             (x, y) float64 32B -99.83 -99.32 -99.79 -99.23
            lat             (x, y) float64 32B 42.25 42.21 42.63 42.59
          * time            (time) datetime64[ns] 24B 2014-09-06 2014-09-07 2014-09-08
            reference_time  datetime64[ns] 8B 2014-09-05
        Dimensions without coordinates: x, y
        Data variables:
            temperature     (x, y, time) float64 96B 29.11 18.2 22.83 ... 16.15 26.63
            precipitation   (x, y, time) float64 96B 5.68 9.256 0.7104 ... 4.615 7.805

        """
        selection = []
        for var_name, variable in self.variables.items():
            has_value_flag = False
            for attr_name, pattern in kwargs.items():
                attr_value = variable.attrs.get(attr_name)
                if (callable(pattern) and pattern(attr_value)) or attr_value == pattern:
                    has_value_flag = True
                else:
                    has_value_flag = False
                    break
            if has_value_flag is True:
                selection.append(var_name)
        return self[selection]

    def unify_chunks(self) -> Self:
        """Unify chunk size along all chunked dimensions of this Dataset.

        Returns
        -------
        Dataset with consistent chunk sizes for all dask-array variables

        See Also
        --------
        dask.array.core.unify_chunks
        """

        return unify_chunks(self)[0]

    def map_blocks(
        self,
        func: Callable[..., T_Xarray],
        args: Sequence[Any] = (),
        kwargs: Mapping[str, Any] | None = None,
        template: DataArray | Dataset | None = None,
    ) -> T_Xarray:
        """
        Apply a function to each block of this Dataset.

        .. warning::
            This method is experimental and its signature may change.

        Parameters
        ----------
        func : callable
            User-provided function that accepts a Dataset as its first
            parameter. The function will receive a subset or 'block' of this Dataset (see below),
            corresponding to one chunk along each chunked dimension. ``func`` will be
            executed as ``func(subset_dataset, *subset_args, **kwargs)``.

            This function must return either a single DataArray or a single Dataset.

            This function cannot add a new chunked dimension.
        args : sequence
            Passed to func after unpacking and subsetting any xarray objects by blocks.
            xarray objects in args must be aligned with obj, otherwise an error is raised.
        kwargs : Mapping or None
            Passed verbatim to func after unpacking. xarray objects, if any, will not be
            subset to blocks. Passing dask collections in kwargs is not allowed.
        template : DataArray, Dataset or None, optional
            xarray object representing the final result after compute is called. If not provided,
            the function will be first run on mocked-up data, that looks like this object but
            has sizes 0, to determine properties of the returned object such as dtype,
            variable names, attributes, new dimensions and new indexes (if any).
            ``template`` must be provided if the function changes the size of existing dimensions.
            When provided, ``attrs`` on variables in `template` are copied over to the result. Any
            ``attrs`` set by ``func`` will be ignored.

        Returns
        -------
        A single DataArray or Dataset with dask backend, reassembled from the outputs of the
        function.

        Notes
        -----
        This function is designed for when ``func`` needs to manipulate a whole xarray object
        subset to each block. Each block is loaded into memory. In the more common case where
        ``func`` can work on numpy arrays, it is recommended to use ``apply_ufunc``.

        If none of the variables in this object is backed by dask arrays, calling this function is
        equivalent to calling ``func(obj, *args, **kwargs)``.

        See Also
        --------
        :func:`dask.array.map_blocks <dask.array.map_blocks>`
        :func:`xarray.apply_ufunc <apply_ufunc>`
        :func:`xarray.DataArray.map_blocks <xarray.DataArray.map_blocks>`

        :doc:`xarray-tutorial:advanced/map_blocks/map_blocks`
            Advanced Tutorial on map_blocks with dask


        Examples
        --------
        Calculate an anomaly from climatology using ``.groupby()``. Using
        ``xr.map_blocks()`` allows for parallel operations with knowledge of ``xarray``,
        its indices, and its methods like ``.groupby()``.

        >>> def calculate_anomaly(da, groupby_type="time.month"):
        ...     gb = da.groupby(groupby_type)
        ...     clim = gb.mean(dim="time")
        ...     return gb - clim
        ...
        >>> time = xr.date_range("1990-01", "1992-01", freq="ME", use_cftime=True)
        >>> month = xr.DataArray(time.month, coords={"time": time}, dims=["time"])
        >>> np.random.seed(123)
        >>> array = xr.DataArray(
        ...     np.random.rand(len(time)),
        ...     dims=["time"],
        ...     coords={"time": time, "month": month},
        ... ).chunk()
        >>> ds = xr.Dataset({"a": array})
        >>> ds.map_blocks(calculate_anomaly, template=ds).compute()
        <xarray.Dataset> Size: 576B
        Dimensions:  (time: 24)
        Coordinates:
          * time     (time) object 192B 1990-01-31 00:00:00 ... 1991-12-31 00:00:00
            month    (time) int64 192B 1 2 3 4 5 6 7 8 9 10 ... 3 4 5 6 7 8 9 10 11 12
        Data variables:
            a        (time) float64 192B 0.1289 0.1132 -0.0856 ... 0.1906 -0.05901

        Note that one must explicitly use ``args=[]`` and ``kwargs={}`` to pass arguments
        to the function being applied in ``xr.map_blocks()``:

        >>> ds.map_blocks(
        ...     calculate_anomaly,
        ...     kwargs={"groupby_type": "time.year"},
        ...     template=ds,
        ... )
        <xarray.Dataset> Size: 576B
        Dimensions:  (time: 24)
        Coordinates:
          * time     (time) object 192B 1990-01-31 00:00:00 ... 1991-12-31 00:00:00
            month    (time) int64 192B dask.array<chunksize=(24,), meta=np.ndarray>
        Data variables:
            a        (time) float64 192B dask.array<chunksize=(24,), meta=np.ndarray>
        """
        from xarray.core.parallel import map_blocks

        return map_blocks(func, self, args, kwargs, template)

    def polyfit(
        self,
        dim: Hashable,
        deg: int,
        skipna: bool | None = None,
        rcond: float | None = None,
        w: Hashable | Any = None,
        full: bool = False,
        cov: bool | Literal["unscaled"] = False,
    ) -> Self:
        """
        Least squares polynomial fit.

        This replicates the behaviour of `numpy.polyfit` but differs by skipping
        invalid values when `skipna = True`.

        Parameters
        ----------
        dim : hashable
            Coordinate along which to fit the polynomials.
        deg : int
            Degree of the fitting polynomial.
        skipna : bool or None, optional
            If True, removes all invalid values before fitting each 1D slices of the array.
            Default is True if data is stored in a dask.array or if there is any
            invalid values, False otherwise.
        rcond : float or None, optional
            Relative condition number to the fit.
        w : hashable or Any, optional
            Weights to apply to the y-coordinate of the sample points.
            Can be an array-like object or the name of a coordinate in the dataset.
        full : bool, default: False
            Whether to return the residuals, matrix rank and singular values in addition
            to the coefficients.
        cov : bool or "unscaled", default: False
            Whether to return to the covariance matrix in addition to the coefficients.
            The matrix is not scaled if `cov='unscaled'`.

        Returns
        -------
        polyfit_results : Dataset
            A single dataset which contains (for each "var" in the input dataset):

            [var]_polyfit_coefficients
                The coefficients of the best fit for each variable in this dataset.
            [var]_polyfit_residuals
                The residuals of the least-square computation for each variable (only included if `full=True`)
                When the matrix rank is deficient, np.nan is returned.
            [dim]_matrix_rank
                The effective rank of the scaled Vandermonde coefficient matrix (only included if `full=True`)
                The rank is computed ignoring the NaN values that might be skipped.
            [dim]_singular_values
                The singular values of the scaled Vandermonde coefficient matrix (only included if `full=True`)
            [var]_polyfit_covariance
                The covariance matrix of the polynomial coefficient estimates (only included if `full=False` and `cov=True`)

        Warns
        -----
        RankWarning
            The rank of the coefficient matrix in the least-squares fit is deficient.
            The warning is not raised with in-memory (not dask) data and `full=True`.

        See Also
        --------
        numpy.polyfit
        numpy.polyval
        xarray.polyval
        """
        from xarray.computation.fit import polyfit as polyfit_impl

        return polyfit_impl(self, dim, deg, skipna, rcond, w, full, cov)

    def pad(
        self,
        pad_width: Mapping[Any, int | tuple[int, int]] | None = None,
        mode: PadModeOptions = "constant",
        stat_length: (
            int | tuple[int, int] | Mapping[Any, tuple[int, int]] | None
        ) = None,
        constant_values: T_DatasetPadConstantValues | None = None,
        end_values: int | tuple[int, int] | Mapping[Any, tuple[int, int]] | None = None,
        reflect_type: PadReflectOptions = None,
        keep_attrs: bool | None = None,
        **pad_width_kwargs: Any,
    ) -> Self:
        """Pad this dataset along one or more dimensions.

        .. warning::
            This function is experimental and its behaviour is likely to change
            especially regarding padding of dimension coordinates (or IndexVariables).

        When using one of the modes ("edge", "reflect", "symmetric", "wrap"),
        coordinates will be padded with the same mode, otherwise coordinates
        are padded using the "constant" mode with fill_value dtypes.NA.

        Parameters
        ----------
        pad_width : mapping of hashable to tuple of int
            Mapping with the form of {dim: (pad_before, pad_after)}
            describing the number of values padded along each dimension.
            {dim: pad} is a shortcut for pad_before = pad_after = pad
        mode : {"constant", "edge", "linear_ramp", "maximum", "mean", "median", \
            "minimum", "reflect", "symmetric", "wrap"}, default: "constant"
            How to pad the DataArray (taken from numpy docs):

            - "constant": Pads with a constant value.
            - "edge": Pads with the edge values of array.
            - "linear_ramp": Pads with the linear ramp between end_value and the
              array edge value.
            - "maximum": Pads with the maximum value of all or part of the
              vector along each axis.
            - "mean": Pads with the mean value of all or part of the
              vector along each axis.
            - "median": Pads with the median value of all or part of the
              vector along each axis.
            - "minimum": Pads with the minimum value of all or part of the
              vector along each axis.
            - "reflect": Pads with the reflection of the vector mirrored on
              the first and last values of the vector along each axis.
            - "symmetric": Pads with the reflection of the vector mirrored
              along the edge of the array.
            - "wrap": Pads with the wrap of the vector along the axis.
              The first values are used to pad the end and the
              end values are used to pad the beginning.

        stat_length : int, tuple or mapping of hashable to tuple, default: None
            Used in 'maximum', 'mean', 'median', and 'minimum'.  Number of
            values at edge of each axis used to calculate the statistic value.
            {dim_1: (before_1, after_1), ... dim_N: (before_N, after_N)} unique
            statistic lengths along each dimension.
            ((before, after),) yields same before and after statistic lengths
            for each dimension.
            (stat_length,) or int is a shortcut for before = after = statistic
            length for all axes.
            Default is ``None``, to use the entire axis.
        constant_values : scalar, tuple, mapping of dim name to scalar or tuple, or \
            mapping of var name to scalar, tuple or to mapping of dim name to scalar or tuple, default: None
            Used in 'constant'. The values to set the padded values for each data variable / axis.
            ``{var_1: {dim_1: (before_1, after_1), ... dim_N: (before_N, after_N)}, ...
            var_M: (before, after)}`` unique pad constants per data variable.
            ``{dim_1: (before_1, after_1), ... dim_N: (before_N, after_N)}`` unique
            pad constants along each dimension.
            ``((before, after),)`` yields same before and after constants for each
            dimension.
            ``(constant,)`` or ``constant`` is a shortcut for ``before = after = constant`` for
            all dimensions.
            Default is ``None``, pads with ``np.nan``.
        end_values : scalar, tuple or mapping of hashable to tuple, default: None
            Used in 'linear_ramp'.  The values used for the ending value of the
            linear_ramp and that will form the edge of the padded array.
            ``{dim_1: (before_1, after_1), ... dim_N: (before_N, after_N)}`` unique
            end values along each dimension.
            ``((before, after),)`` yields same before and after end values for each
            axis.
            ``(constant,)`` or ``constant`` is a shortcut for ``before = after = constant`` for
            all axes.
            Default is None.
        reflect_type : {"even", "odd", None}, optional
            Used in "reflect", and "symmetric".  The "even" style is the
            default with an unaltered reflection around the edge value.  For
            the "odd" style, the extended part of the array is created by
            subtracting the reflected values from two times the edge value.
        keep_attrs : bool or None, optional
            If True, the attributes (``attrs``) will be copied from the
            original object to the new one. If False, the new object
            will be returned without attributes.
        **pad_width_kwargs
            The keyword arguments form of ``pad_width``.
            One of ``pad_width`` or ``pad_width_kwargs`` must be provided.

        Returns
        -------
        padded : Dataset
            Dataset with the padded coordinates and data.

        See Also
        --------
        Dataset.shift, Dataset.roll, Dataset.bfill, Dataset.ffill, numpy.pad, dask.array.pad

        Notes
        -----
        By default when ``mode="constant"`` and ``constant_values=None``, integer types will be
        promoted to ``float`` and padded with ``np.nan``. To avoid type promotion
        specify ``constant_values=np.nan``

        Padding coordinates will drop their corresponding index (if any) and will reset default
        indexes for dimension coordinates.

        Examples
        --------
        >>> ds = xr.Dataset({"foo": ("x", range(5))})
        >>> ds.pad(x=(1, 2))
        <xarray.Dataset> Size: 64B
        Dimensions:  (x: 8)
        Dimensions without coordinates: x
        Data variables:
            foo      (x) float64 64B nan 0.0 1.0 2.0 3.0 4.0 nan nan
        """
        pad_width = either_dict_or_kwargs(pad_width, pad_width_kwargs, "pad")

        if mode in ("edge", "reflect", "symmetric", "wrap"):
            coord_pad_mode = mode
            coord_pad_options = {
                "stat_length": stat_length,
                "constant_values": constant_values,
                "end_values": end_values,
                "reflect_type": reflect_type,
            }
        else:
            coord_pad_mode = "constant"
            coord_pad_options = {}

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)

        variables = {}

        # keep indexes that won't be affected by pad and drop all other indexes
        xindexes = self.xindexes
        pad_dims = set(pad_width)
        indexes = {
            k: idx
            for k, idx in xindexes.items()
            if not pad_dims.intersection(xindexes.get_all_dims(k))
        }

        for name, var in self.variables.items():
            var_pad_width = {k: v for k, v in pad_width.items() if k in var.dims}
            if not var_pad_width:
                variables[name] = var
            elif name in self.data_vars:
                if utils.is_dict_like(constant_values):
                    if name in constant_values.keys():
                        filtered_constant_values = constant_values[name]
                    elif not set(var.dims).isdisjoint(constant_values.keys()):
                        filtered_constant_values = {
                            k: v for k, v in constant_values.items() if k in var.dims
                        }
                    else:
                        filtered_constant_values = 0  # TODO: https://github.com/pydata/xarray/pull/9353#discussion_r1724018352
                else:
                    filtered_constant_values = constant_values
                variables[name] = var.pad(
                    pad_width=var_pad_width,
                    mode=mode,
                    stat_length=stat_length,
                    constant_values=filtered_constant_values,
                    end_values=end_values,
                    reflect_type=reflect_type,
                    keep_attrs=keep_attrs,
                )
            else:
                variables[name] = var.pad(
                    pad_width=var_pad_width,
                    mode=coord_pad_mode,
                    keep_attrs=keep_attrs,
                    **coord_pad_options,  # type: ignore[arg-type]
                )
                # reset default index of dimension coordinates
                if (name,) == var.dims:
                    dim_var = {name: variables[name]}
                    index = PandasIndex.from_variables(dim_var, options={})
                    index_vars = index.create_variables(dim_var)
                    indexes[name] = index
                    variables[name] = index_vars[name]

        attrs = self._attrs if keep_attrs else None
        return self._replace_with_new_dims(variables, indexes=indexes, attrs=attrs)

    def idxmin(
        self,
        dim: Hashable | None = None,
        *,
        skipna: bool | None = None,
        fill_value: Any = xrdtypes.NA,
        keep_attrs: bool | None = None,
    ) -> Self:
        """Return the coordinate label of the minimum value along a dimension.

        Returns a new `Dataset` named after the dimension with the values of
        the coordinate labels along that dimension corresponding to minimum
        values along that dimension.

        In comparison to :py:meth:`~Dataset.argmin`, this returns the
        coordinate label while :py:meth:`~Dataset.argmin` returns the index.

        Parameters
        ----------
        dim : Hashable, optional
            Dimension over which to apply `idxmin`.  This is optional for 1D
            variables, but required for variables with 2 or more dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for ``float``, ``complex``, and ``object``
            dtypes; other dtypes either do not have a sentinel missing value
            (``int``) or ``skipna=True`` has not been implemented
            (``datetime64`` or ``timedelta64``).
        fill_value : Any, default: NaN
            Value to be filled in case all of the values along a dimension are
            null.  By default this is NaN.  The fill value and result are
            automatically converted to a compatible dtype if possible.
            Ignored if ``skipna`` is False.
        keep_attrs : bool or None, optional
            If True, the attributes (``attrs``) will be copied from the
            original object to the new one. If False, the new object
            will be returned without attributes.

        Returns
        -------
        reduced : Dataset
            New `Dataset` object with `idxmin` applied to its data and the
            indicated dimension removed.

        See Also
        --------
        DataArray.idxmin, Dataset.idxmax, Dataset.min, Dataset.argmin

        Examples
        --------
        >>> array1 = xr.DataArray(
        ...     [0, 2, 1, 0, -2], dims="x", coords={"x": ["a", "b", "c", "d", "e"]}
        ... )
        >>> array2 = xr.DataArray(
        ...     [
        ...         [2.0, 1.0, 2.0, 0.0, -2.0],
        ...         [-4.0, np.nan, 2.0, np.nan, -2.0],
        ...         [np.nan, np.nan, 1.0, np.nan, np.nan],
        ...     ],
        ...     dims=["y", "x"],
        ...     coords={"y": [-1, 0, 1], "x": ["a", "b", "c", "d", "e"]},
        ... )
        >>> ds = xr.Dataset({"int": array1, "float": array2})
        >>> ds.min(dim="x")
        <xarray.Dataset> Size: 56B
        Dimensions:  (y: 3)
        Coordinates:
          * y        (y) int64 24B -1 0 1
        Data variables:
            int      int64 8B -2
            float    (y) float64 24B -2.0 -4.0 1.0
        >>> ds.argmin(dim="x")
        <xarray.Dataset> Size: 56B
        Dimensions:  (y: 3)
        Coordinates:
          * y        (y) int64 24B -1 0 1
        Data variables:
            int      int64 8B 4
            float    (y) int64 24B 4 0 2
        >>> ds.idxmin(dim="x")
        <xarray.Dataset> Size: 52B
        Dimensions:  (y: 3)
        Coordinates:
          * y        (y) int64 24B -1 0 1
        Data variables:
            int      <U1 4B 'e'
            float    (y) object 24B 'e' 'a' 'c'
        """
        return self.map(
            methodcaller(
                "idxmin",
                dim=dim,
                skipna=skipna,
                fill_value=fill_value,
                keep_attrs=keep_attrs,
            )
        )

    def idxmax(
        self,
        dim: Hashable | None = None,
        *,
        skipna: bool | None = None,
        fill_value: Any = xrdtypes.NA,
        keep_attrs: bool | None = None,
    ) -> Self:
        """Return the coordinate label of the maximum value along a dimension.

        Returns a new `Dataset` named after the dimension with the values of
        the coordinate labels along that dimension corresponding to maximum
        values along that dimension.

        In comparison to :py:meth:`~Dataset.argmax`, this returns the
        coordinate label while :py:meth:`~Dataset.argmax` returns the index.

        Parameters
        ----------
        dim : str, optional
            Dimension over which to apply `idxmax`.  This is optional for 1D
            variables, but required for variables with 2 or more dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for ``float``, ``complex``, and ``object``
            dtypes; other dtypes either do not have a sentinel missing value
            (``int``) or ``skipna=True`` has not been implemented
            (``datetime64`` or ``timedelta64``).
        fill_value : Any, default: NaN
            Value to be filled in case all of the values along a dimension are
            null.  By default this is NaN.  The fill value and result are
            automatically converted to a compatible dtype if possible.
            Ignored if ``skipna`` is False.
        keep_attrs : bool or None, optional
            If True, the attributes (``attrs``) will be copied from the
            original object to the new one. If False, the new object
            will be returned without attributes.

        Returns
        -------
        reduced : Dataset
            New `Dataset` object with `idxmax` applied to its data and the
            indicated dimension removed.

        See Also
        --------
        DataArray.idxmax, Dataset.idxmin, Dataset.max, Dataset.argmax

        Examples
        --------
        >>> array1 = xr.DataArray(
        ...     [0, 2, 1, 0, -2], dims="x", coords={"x": ["a", "b", "c", "d", "e"]}
        ... )
        >>> array2 = xr.DataArray(
        ...     [
        ...         [2.0, 1.0, 2.0, 0.0, -2.0],
        ...         [-4.0, np.nan, 2.0, np.nan, -2.0],
        ...         [np.nan, np.nan, 1.0, np.nan, np.nan],
        ...     ],
        ...     dims=["y", "x"],
        ...     coords={"y": [-1, 0, 1], "x": ["a", "b", "c", "d", "e"]},
        ... )
        >>> ds = xr.Dataset({"int": array1, "float": array2})
        >>> ds.max(dim="x")
        <xarray.Dataset> Size: 56B
        Dimensions:  (y: 3)
        Coordinates:
          * y        (y) int64 24B -1 0 1
        Data variables:
            int      int64 8B 2
            float    (y) float64 24B 2.0 2.0 1.0
        >>> ds.argmax(dim="x")
        <xarray.Dataset> Size: 56B
        Dimensions:  (y: 3)
        Coordinates:
          * y        (y) int64 24B -1 0 1
        Data variables:
            int      int64 8B 1
            float    (y) int64 24B 0 2 2
        >>> ds.idxmax(dim="x")
        <xarray.Dataset> Size: 52B
        Dimensions:  (y: 3)
        Coordinates:
          * y        (y) int64 24B -1 0 1
        Data variables:
            int      <U1 4B 'b'
            float    (y) object 24B 'a' 'c' 'c'
        """
        return self.map(
            methodcaller(
                "idxmax",
                dim=dim,
                skipna=skipna,
                fill_value=fill_value,
                keep_attrs=keep_attrs,
            )
        )

    def argmin(self, dim: Hashable | None = None, **kwargs) -> Self:
        """Indices of the minima of the member variables.

        If there are multiple minima, the indices of the first one found will be
        returned.

        Parameters
        ----------
        dim : Hashable, optional
            The dimension over which to find the minimum. By default, finds minimum over
            all dimensions - for now returning an int for backward compatibility, but
            this is deprecated, in future will be an error, since DataArray.argmin will
            return a dict with indices for all dimensions, which does not make sense for
            a Dataset.
        keep_attrs : bool, optional
            If True, the attributes (`attrs`) will be copied from the original
            object to the new one.  If False (default), the new object will be
            returned without attributes.
        skipna : bool, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or skipna=True has not been
            implemented (object, datetime64 or timedelta64).

        Returns
        -------
        result : Dataset

        Examples
        --------
        >>> dataset = xr.Dataset(
        ...     {
        ...         "math_scores": (
        ...             ["student", "test"],
        ...             [[90, 85, 79], [78, 80, 85], [95, 92, 98]],
        ...         ),
        ...         "english_scores": (
        ...             ["student", "test"],
        ...             [[88, 90, 92], [75, 82, 79], [39, 96, 78]],
        ...         ),
        ...     },
        ...     coords={
        ...         "student": ["Alice", "Bob", "Charlie"],
        ...         "test": ["Test 1", "Test 2", "Test 3"],
        ...     },
        ... )

        # Indices of the minimum values along the 'student' dimension are calculated

        >>> argmin_indices = dataset.argmin(dim="student")

        >>> min_score_in_math = dataset["student"].isel(
        ...     student=argmin_indices["math_scores"]
        ... )
        >>> min_score_in_math
        <xarray.DataArray 'student' (test: 3)> Size: 84B
        array(['Bob', 'Bob', 'Alice'], dtype='<U7')
        Coordinates:
            student  (test) <U7 84B 'Bob' 'Bob' 'Alice'
          * test     (test) <U6 72B 'Test 1' 'Test 2' 'Test 3'

        >>> min_score_in_english = dataset["student"].isel(
        ...     student=argmin_indices["english_scores"]
        ... )
        >>> min_score_in_english
        <xarray.DataArray 'student' (test: 3)> Size: 84B
        array(['Charlie', 'Bob', 'Charlie'], dtype='<U7')
        Coordinates:
            student  (test) <U7 84B 'Charlie' 'Bob' 'Charlie'
          * test     (test) <U6 72B 'Test 1' 'Test 2' 'Test 3'

        See Also
        --------
        Dataset.idxmin
        DataArray.argmin
        """
        if dim is None:
            warnings.warn(
                "Once the behaviour of DataArray.argmin() and Variable.argmin() without "
                "dim changes to return a dict of indices of each dimension, for "
                "consistency it will be an error to call Dataset.argmin() with no argument,"
                "since we don't return a dict of Datasets.",
                DeprecationWarning,
                stacklevel=2,
            )
        if (
            dim is None
            or (not isinstance(dim, Sequence) and dim is not ...)
            or isinstance(dim, str)
        ):
            # Return int index if single dimension is passed, and is not part of a
            # sequence
            argmin_func = duck_array_ops.argmin
            return self.reduce(
                argmin_func, dim=None if dim is None else [dim], **kwargs
            )
        else:
            raise ValueError(
                "When dim is a sequence or ..., DataArray.argmin() returns a dict. "
                "dicts cannot be contained in a Dataset, so cannot call "
                "Dataset.argmin() with a sequence or ... for dim"
            )

    def argmax(self, dim: Hashable | None = None, **kwargs) -> Self:
        """Indices of the maxima of the member variables.

        If there are multiple maxima, the indices of the first one found will be
        returned.

        Parameters
        ----------
        dim : str, optional
            The dimension over which to find the maximum. By default, finds maximum over
            all dimensions - for now returning an int for backward compatibility, but
            this is deprecated, in future will be an error, since DataArray.argmax will
            return a dict with indices for all dimensions, which does not make sense for
            a Dataset.
        keep_attrs : bool, optional
            If True, the attributes (`attrs`) will be copied from the original
            object to the new one.  If False (default), the new object will be
            returned without attributes.
        skipna : bool, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or skipna=True has not been
            implemented (object, datetime64 or timedelta64).

        Returns
        -------
        result : Dataset

        Examples
        --------

        >>> dataset = xr.Dataset(
        ...     {
        ...         "math_scores": (
        ...             ["student", "test"],
        ...             [[90, 85, 92], [78, 80, 85], [95, 92, 98]],
        ...         ),
        ...         "english_scores": (
        ...             ["student", "test"],
        ...             [[88, 90, 92], [75, 82, 79], [93, 96, 91]],
        ...         ),
        ...     },
        ...     coords={
        ...         "student": ["Alice", "Bob", "Charlie"],
        ...         "test": ["Test 1", "Test 2", "Test 3"],
        ...     },
        ... )

        # Indices of the maximum values along the 'student' dimension are calculated

        >>> argmax_indices = dataset.argmax(dim="test")

        >>> argmax_indices
        <xarray.Dataset> Size: 132B
        Dimensions:         (student: 3)
        Coordinates:
          * student         (student) <U7 84B 'Alice' 'Bob' 'Charlie'
        Data variables:
            math_scores     (student) int64 24B 2 2 2
            english_scores  (student) int64 24B 2 1 1

        See Also
        --------
        DataArray.argmax

        """
        if dim is None:
            warnings.warn(
                "Once the behaviour of DataArray.argmin() and Variable.argmin() without "
                "dim changes to return a dict of indices of each dimension, for "
                "consistency it will be an error to call Dataset.argmin() with no argument,"
                "since we don't return a dict of Datasets.",
                DeprecationWarning,
                stacklevel=2,
            )
        if (
            dim is None
            or (not isinstance(dim, Sequence) and dim is not ...)
            or isinstance(dim, str)
        ):
            # Return int index if single dimension is passed, and is not part of a
            # sequence
            argmax_func = duck_array_ops.argmax
            return self.reduce(
                argmax_func, dim=None if dim is None else [dim], **kwargs
            )
        else:
            raise ValueError(
                "When dim is a sequence or ..., DataArray.argmin() returns a dict. "
                "dicts cannot be contained in a Dataset, so cannot call "
                "Dataset.argmin() with a sequence or ... for dim"
            )

    def eval(
        self,
        statement: str,
        *,
        parser: QueryParserOptions = "pandas",
    ) -> Self | T_DataArray:
        """
        Calculate an expression supplied as a string in the context of the dataset.

        This is currently experimental; the API may change particularly around
        assignments, which currently return a ``Dataset`` with the additional variable.
        Currently only the ``python`` engine is supported, which has the same
        performance as executing in python.

        Parameters
        ----------
        statement : str
            String containing the Python-like expression to evaluate.

        Returns
        -------
        result : Dataset or DataArray, depending on whether ``statement`` contains an
        assignment.

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     {"a": ("x", np.arange(0, 5, 1)), "b": ("x", np.linspace(0, 1, 5))}
        ... )
        >>> ds
        <xarray.Dataset> Size: 80B
        Dimensions:  (x: 5)
        Dimensions without coordinates: x
        Data variables:
            a        (x) int64 40B 0 1 2 3 4
            b        (x) float64 40B 0.0 0.25 0.5 0.75 1.0

        >>> ds.eval("a + b")
        <xarray.DataArray (x: 5)> Size: 40B
        array([0.  , 1.25, 2.5 , 3.75, 5.  ])
        Dimensions without coordinates: x

        >>> ds.eval("c = a + b")
        <xarray.Dataset> Size: 120B
        Dimensions:  (x: 5)
        Dimensions without coordinates: x
        Data variables:
            a        (x) int64 40B 0 1 2 3 4
            b        (x) float64 40B 0.0 0.25 0.5 0.75 1.0
            c        (x) float64 40B 0.0 1.25 2.5 3.75 5.0
        """

        return pd.eval(  # type: ignore[return-value]
            statement,
            resolvers=[self],
            target=self,
            parser=parser,
            # Because numexpr returns a numpy array, using that engine results in
            # different behavior. We'd be very open to a contribution handling this.
            engine="python",
        )

    def query(
        self,
        queries: Mapping[Any, Any] | None = None,
        parser: QueryParserOptions = "pandas",
        engine: QueryEngineOptions = None,
        missing_dims: ErrorOptionsWithWarn = "raise",
        **queries_kwargs: Any,
    ) -> Self:
        """Return a new dataset with each array indexed along the specified
        dimension(s), where the indexers are given as strings containing
        Python expressions to be evaluated against the data variables in the
        dataset.

        Parameters
        ----------
        queries : dict-like, optional
            A dict-like with keys matching dimensions and values given by strings
            containing Python expressions to be evaluated against the data variables
            in the dataset. The expressions will be evaluated using the pandas
            eval() function, and can contain any valid Python expressions but cannot
            contain any Python statements.
        parser : {"pandas", "python"}, default: "pandas"
            The parser to use to construct the syntax tree from the expression.
            The default of 'pandas' parses code slightly different than standard
            Python. Alternatively, you can parse an expression using the 'python'
            parser to retain strict Python semantics.
        engine : {"python", "numexpr", None}, default: None
            The engine used to evaluate the expression. Supported engines are:

            - None: tries to use numexpr, falls back to python
            - "numexpr": evaluates expressions using numexpr
            - "python": performs operations as if you had eval’d in top level python

        missing_dims : {"raise", "warn", "ignore"}, default: "raise"
            What to do if dimensions that should be selected from are not present in the
            Dataset:

            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        **queries_kwargs : {dim: query, ...}, optional
            The keyword arguments form of ``queries``.
            One of queries or queries_kwargs must be provided.

        Returns
        -------
        obj : Dataset
            A new Dataset with the same contents as this dataset, except each
            array and dimension is indexed by the results of the appropriate
            queries.

        See Also
        --------
        Dataset.isel
        pandas.eval

        Examples
        --------
        >>> a = np.arange(0, 5, 1)
        >>> b = np.linspace(0, 1, 5)
        >>> ds = xr.Dataset({"a": ("x", a), "b": ("x", b)})
        >>> ds
        <xarray.Dataset> Size: 80B
        Dimensions:  (x: 5)
        Dimensions without coordinates: x
        Data variables:
            a        (x) int64 40B 0 1 2 3 4
            b        (x) float64 40B 0.0 0.25 0.5 0.75 1.0
        >>> ds.query(x="a > 2")
        <xarray.Dataset> Size: 32B
        Dimensions:  (x: 2)
        Dimensions without coordinates: x
        Data variables:
            a        (x) int64 16B 3 4
            b        (x) float64 16B 0.75 1.0
        """

        # allow queries to be given either as a dict or as kwargs
        queries = either_dict_or_kwargs(queries, queries_kwargs, "query")

        # check queries
        for dim, expr in queries.items():
            if not isinstance(expr, str):
                msg = f"expr for dim {dim} must be a string to be evaluated, {type(expr)} given"
                raise ValueError(msg)

        # evaluate the queries to create the indexers
        indexers = {
            dim: pd.eval(expr, resolvers=[self], parser=parser, engine=engine)
            for dim, expr in queries.items()
        }

        # apply the selection
        return self.isel(indexers, missing_dims=missing_dims)

    def curvefit(
        self,
        coords: str | DataArray | Iterable[str | DataArray],
        func: Callable[..., Any],
        reduce_dims: Dims = None,
        skipna: bool = True,
        p0: Mapping[str, float | DataArray] | None = None,
        bounds: Mapping[str, tuple[float | DataArray, float | DataArray]] | None = None,
        param_names: Sequence[str] | None = None,
        errors: ErrorOptions = "raise",
        kwargs: dict[str, Any] | None = None,
    ) -> Self:
        """
        Curve fitting optimization for arbitrary functions.

        Wraps :py:func:`scipy.optimize.curve_fit` with :py:func:`~xarray.apply_ufunc`.

        Parameters
        ----------
        coords : hashable, DataArray, or sequence of hashable or DataArray
            Independent coordinate(s) over which to perform the curve fitting. Must share
            at least one dimension with the calling object. When fitting multi-dimensional
            functions, supply `coords` as a sequence in the same order as arguments in
            `func`. To fit along existing dimensions of the calling object, `coords` can
            also be specified as a str or sequence of strs.
        func : callable
            User specified function in the form `f(x, *params)` which returns a numpy
            array of length `len(x)`. `params` are the fittable parameters which are optimized
            by scipy curve_fit. `x` can also be specified as a sequence containing multiple
            coordinates, e.g. `f((x0, x1), *params)`.
        reduce_dims : str, Iterable of Hashable or None, optional
            Additional dimension(s) over which to aggregate while fitting. For example,
            calling `ds.curvefit(coords='time', reduce_dims=['lat', 'lon'], ...)` will
            aggregate all lat and lon points and fit the specified function along the
            time dimension.
        skipna : bool, default: True
            Whether to skip missing values when fitting. Default is True.
        p0 : dict-like, optional
            Optional dictionary of parameter names to initial guesses passed to the
            `curve_fit` `p0` arg. If the values are DataArrays, they will be appropriately
            broadcast to the coordinates of the array. If none or only some parameters are
            passed, the rest will be assigned initial values following the default scipy
            behavior.
        bounds : dict-like, optional
            Optional dictionary of parameter names to tuples of bounding values passed to the
            `curve_fit` `bounds` arg. If any of the bounds are DataArrays, they will be
            appropriately broadcast to the coordinates of the array. If none or only some
            parameters are passed, the rest will be unbounded following the default scipy
            behavior.
        param_names : sequence of hashable, optional
            Sequence of names for the fittable parameters of `func`. If not supplied,
            this will be automatically determined by arguments of `func`. `param_names`
            should be manually supplied when fitting a function that takes a variable
            number of parameters.
        errors : {"raise", "ignore"}, default: "raise"
            If 'raise', any errors from the `scipy.optimize_curve_fit` optimization will
            raise an exception. If 'ignore', the coefficients and covariances for the
            coordinates where the fitting failed will be NaN.
        **kwargs : optional
            Additional keyword arguments to passed to scipy curve_fit.

        Returns
        -------
        curvefit_results : Dataset
            A single dataset which contains:

            [var]_curvefit_coefficients
                The coefficients of the best fit.
            [var]_curvefit_covariance
                The covariance matrix of the coefficient estimates.

        See Also
        --------
        Dataset.polyfit
        scipy.optimize.curve_fit
        xarray.Dataset.xlm.modelfit
            External method from `xarray-lmfit <https://xarray-lmfit.readthedocs.io/>`_
            with more curve fitting functionality.
        """
        from xarray.computation.fit import curvefit as curvefit_impl

        return curvefit_impl(
            self,
            coords,
            func,
            reduce_dims,
            skipna,
            p0,
            bounds,
            param_names,
            errors,
            kwargs,
        )

    def drop_duplicates(
        self,
        dim: Hashable | Iterable[Hashable],
        *,
        keep: Literal["first", "last", False] = "first",
    ) -> Self:
        """Returns a new Dataset with duplicate dimension values removed.

        Parameters
        ----------
        dim : dimension label or labels
            Pass `...` to drop duplicates along all dimensions.
        keep : {"first", "last", False}, default: "first"
            Determines which duplicates (if any) to keep.
            - ``"first"`` : Drop duplicates except for the first occurrence.
            - ``"last"`` : Drop duplicates except for the last occurrence.
            - False : Drop all duplicates.

        Returns
        -------
        Dataset

        See Also
        --------
        DataArray.drop_duplicates
        """
        if isinstance(dim, str):
            dims: Iterable = (dim,)
        elif dim is ...:
            dims = self.dims
        elif not isinstance(dim, Iterable):
            dims = [dim]
        else:
            dims = dim

        missing_dims = set(dims) - set(self.dims)
        if missing_dims:
            raise ValueError(
                f"Dimensions {tuple(missing_dims)} not found in data dimensions {tuple(self.dims)}"
            )

        indexes = {dim: ~self.get_index(dim).duplicated(keep=keep) for dim in dims}
        return self.isel(indexes)

    def convert_calendar(
        self,
        calendar: CFCalendar,
        dim: Hashable = "time",
        align_on: Literal["date", "year"] | None = None,
        missing: Any | None = None,
        use_cftime: bool | None = None,
    ) -> Self:
        """Convert the Dataset to another calendar.

        Only converts the individual timestamps, does not modify any data except
        in dropping invalid/surplus dates or inserting missing dates.

        If the source and target calendars are either no_leap, all_leap or a
        standard type, only the type of the time array is modified.
        When converting to a leap year from a non-leap year, the 29th of February
        is removed from the array. In the other direction the 29th of February
        will be missing in the output, unless `missing` is specified,
        in which case that value is inserted.

        For conversions involving `360_day` calendars, see Notes.

        This method is safe to use with sub-daily data as it doesn't touch the
        time part of the timestamps.

        Parameters
        ---------
        calendar : str
            The target calendar name.
        dim : Hashable, default: "time"
            Name of the time coordinate.
        align_on : {None, 'date', 'year'}, optional
            Must be specified when either source or target is a `360_day` calendar,
            ignored otherwise. See Notes.
        missing : Any or None, optional
            By default, i.e. if the value is None, this method will simply attempt
            to convert the dates in the source calendar to the same dates in the
            target calendar, and drop any of those that are not possible to
            represent.  If a value is provided, a new time coordinate will be
            created in the target calendar with the same frequency as the original
            time coordinate; for any dates that are not present in the source, the
            data will be filled with this value.  Note that using this mode requires
            that the source data have an inferable frequency; for more information
            see :py:func:`xarray.infer_freq`.  For certain frequency, source, and
            target calendar combinations, this could result in many missing values, see notes.
        use_cftime : bool or None, optional
            Whether to use cftime objects in the output, only used if `calendar`
            is one of {"proleptic_gregorian", "gregorian" or "standard"}.
            If True, the new time axis uses cftime objects.
            If None (default), it uses :py:class:`numpy.datetime64` values if the
            date range permits it, and :py:class:`cftime.datetime` objects if not.
            If False, it uses :py:class:`numpy.datetime64`  or fails.

        Returns
        -------
        Dataset
            Copy of the dataarray with the time coordinate converted to the
            target calendar. If 'missing' was None (default), invalid dates in
            the new calendar are dropped, but missing dates are not inserted.
            If `missing` was given, the new data is reindexed to have a time axis
            with the same frequency as the source, but in the new calendar; any
            missing datapoints are filled with `missing`.

        Notes
        -----
        Passing a value to `missing` is only usable if the source's time coordinate as an
        inferable frequencies (see :py:func:`~xarray.infer_freq`) and is only appropriate
        if the target coordinate, generated from this frequency, has dates equivalent to the
        source. It is usually **not** appropriate to use this mode with:

        - Period-end frequencies : 'A', 'Y', 'Q' or 'M', in opposition to 'AS' 'YS', 'QS' and 'MS'
        - Sub-monthly frequencies that do not divide a day evenly : 'W', 'nD' where `N != 1`
            or 'mH' where 24 % m != 0).

        If one of the source or target calendars is `"360_day"`, `align_on` must
        be specified and two options are offered.

        - "year"
            The dates are translated according to their relative position in the year,
            ignoring their original month and day information, meaning that the
            missing/surplus days are added/removed at regular intervals.

            From a `360_day` to a standard calendar, the output will be missing the
            following dates (day of year in parentheses):

            To a leap year:
                January 31st (31), March 31st (91), June 1st (153), July 31st (213),
                September 31st (275) and November 30th (335).
            To a non-leap year:
                February 6th (36), April 19th (109), July 2nd (183),
                September 12th (255), November 25th (329).

            From a standard calendar to a `"360_day"`, the following dates in the
            source array will be dropped:

            From a leap year:
                January 31st (31), April 1st (92), June 1st (153), August 1st (214),
                September 31st (275), December 1st (336)
            From a non-leap year:
                February 6th (37), April 20th (110), July 2nd (183),
                September 13th (256), November 25th (329)

            This option is best used on daily and subdaily data.

        - "date"
            The month/day information is conserved and invalid dates are dropped
            from the output. This means that when converting from a `"360_day"` to a
            standard calendar, all 31st (Jan, March, May, July, August, October and
            December) will be missing as there is no equivalent dates in the
            `"360_day"` calendar and the 29th (on non-leap years) and 30th of February
            will be dropped as there are no equivalent dates in a standard calendar.

            This option is best used with data on a frequency coarser than daily.
        """
        return convert_calendar(
            self,
            calendar,
            dim=dim,
            align_on=align_on,
            missing=missing,
            use_cftime=use_cftime,
        )

    def interp_calendar(
        self,
        target: pd.DatetimeIndex | CFTimeIndex | DataArray,
        dim: Hashable = "time",
    ) -> Self:
        """Interpolates the Dataset to another calendar based on decimal year measure.

        Each timestamp in `source` and `target` are first converted to their decimal
        year equivalent then `source` is interpolated on the target coordinate.
        The decimal year of a timestamp is its year plus its sub-year component
        converted to the fraction of its year. For example "2000-03-01 12:00" is
        2000.1653 in a standard calendar or 2000.16301 in a `"noleap"` calendar.

        This method should only be used when the time (HH:MM:SS) information of
        time coordinate is not important.

        Parameters
        ----------
        target: DataArray or DatetimeIndex or CFTimeIndex
            The target time coordinate of a valid dtype
            (np.datetime64 or cftime objects)
        dim : Hashable, default: "time"
            The time coordinate name.

        Return
        ------
        DataArray
            The source interpolated on the decimal years of target,
        """
        return interp_calendar(self, target, dim=dim)

    @_deprecate_positional_args("v2024.07.0")
    def groupby(
        self,
        group: GroupInput = None,
        *,
        squeeze: Literal[False] = False,
        restore_coord_dims: bool = False,
        eagerly_compute_group: Literal[False] | None = None,
        **groupers: Grouper,
    ) -> DatasetGroupBy:
        """Returns a DatasetGroupBy object for performing grouped operations.

        Parameters
        ----------
        group : str or DataArray or IndexVariable or sequence of hashable or mapping of hashable to Grouper
            Array whose unique values should be used to group this array. If a
            Hashable, must be the name of a coordinate contained in this dataarray. If a dictionary,
            must map an existing variable name to a :py:class:`Grouper` instance.
        squeeze : False
            This argument is deprecated.
        restore_coord_dims : bool, default: False
            If True, also restore the dimension order of multi-dimensional
            coordinates.
        eagerly_compute_group: False, optional
            This argument is deprecated.
        **groupers : Mapping of str to Grouper or Resampler
            Mapping of variable name to group by to :py:class:`Grouper` or :py:class:`Resampler` object.
            One of ``group`` or ``groupers`` must be provided.
            Only a single ``grouper`` is allowed at present.

        Returns
        -------
        grouped : DatasetGroupBy
            A `DatasetGroupBy` object patterned after `pandas.GroupBy` that can be
            iterated over in the form of `(unique_value, grouped_array)` pairs.

        Examples
        --------
        >>> ds = xr.Dataset(
        ...     {"foo": (("x", "y"), np.arange(12).reshape((4, 3)))},
        ...     coords={"x": [10, 20, 30, 40], "letters": ("x", list("abba"))},
        ... )

        Grouping by a single variable is easy

        >>> ds.groupby("letters")
        <DatasetGroupBy, grouped over 1 grouper(s), 2 groups in total:
            'letters': UniqueGrouper('letters'), 2/2 groups with labels 'a', 'b'>

        Execute a reduction

        >>> ds.groupby("letters").sum()
        <xarray.Dataset> Size: 64B
        Dimensions:  (letters: 2, y: 3)
        Coordinates:
          * letters  (letters) object 16B 'a' 'b'
        Dimensions without coordinates: y
        Data variables:
            foo      (letters, y) int64 48B 9 11 13 9 11 13

        Grouping by multiple variables

        >>> ds.groupby(["letters", "x"])
        <DatasetGroupBy, grouped over 2 grouper(s), 8 groups in total:
            'letters': UniqueGrouper('letters'), 2/2 groups with labels 'a', 'b'
            'x': UniqueGrouper('x'), 4/4 groups with labels 10, 20, 30, 40>

        Use Grouper objects to express more complicated GroupBy operations

        >>> from xarray.groupers import BinGrouper, UniqueGrouper
        >>>
        >>> ds.groupby(x=BinGrouper(bins=[5, 15, 25]), letters=UniqueGrouper()).sum()
        <xarray.Dataset> Size: 144B
        Dimensions:  (y: 3, x_bins: 2, letters: 2)
        Coordinates:
          * x_bins   (x_bins) interval[int64, right] 32B (5, 15] (15, 25]
          * letters  (letters) object 16B 'a' 'b'
        Dimensions without coordinates: y
        Data variables:
            foo      (y, x_bins, letters) float64 96B 0.0 nan nan 3.0 ... nan nan 5.0

        See Also
        --------
        :ref:`groupby`
            Users guide explanation of how to group and bin data.

        :doc:`xarray-tutorial:intermediate/computation/01-high-level-computation-patterns`
            Tutorial on :py:func:`~xarray.Dataset.Groupby` for windowed computation.

        :doc:`xarray-tutorial:fundamentals/03.2_groupby_with_xarray`
            Tutorial on :py:func:`~xarray.Dataset.Groupby` demonstrating reductions, transformation and comparison with :py:func:`~xarray.Dataset.resample`.

        :external:py:meth:`pandas.DataFrame.groupby <pandas.DataFrame.groupby>`
        :func:`Dataset.groupby_bins <Dataset.groupby_bins>`
        :func:`DataArray.groupby <DataArray.groupby>`
        :class:`core.groupby.DatasetGroupBy`
        :func:`Dataset.coarsen <Dataset.coarsen>`
        :func:`Dataset.resample <Dataset.resample>`
        :func:`DataArray.resample <DataArray.resample>`
        """
        from xarray.core.groupby import (
            DatasetGroupBy,
            _parse_group_and_groupers,
            _validate_groupby_squeeze,
        )

        _validate_groupby_squeeze(squeeze)
        rgroupers = _parse_group_and_groupers(
            self, group, groupers, eagerly_compute_group=eagerly_compute_group
        )

        return DatasetGroupBy(self, rgroupers, restore_coord_dims=restore_coord_dims)

    @_deprecate_positional_args("v2024.07.0")
    def groupby_bins(
        self,
        group: Hashable | DataArray | IndexVariable,
        bins: Bins,
        right: bool = True,
        labels: ArrayLike | None = None,
        precision: int = 3,
        include_lowest: bool = False,
        squeeze: Literal[False] = False,
        restore_coord_dims: bool = False,
        duplicates: Literal["raise", "drop"] = "raise",
        eagerly_compute_group: Literal[False] | None = None,
    ) -> DatasetGroupBy:
        """Returns a DatasetGroupBy object for performing grouped operations.

        Rather than using all unique values of `group`, the values are discretized
        first by applying `pandas.cut` [1]_ to `group`.

        Parameters
        ----------
        group : Hashable, DataArray or IndexVariable
            Array whose binned values should be used to group this array. If a
            string, must be the name of a variable contained in this dataset.
        bins : int or array-like
            If bins is an int, it defines the number of equal-width bins in the
            range of x. However, in this case, the range of x is extended by .1%
            on each side to include the min or max values of x. If bins is a
            sequence it defines the bin edges allowing for non-uniform bin
            width. No extension of the range of x is done in this case.
        right : bool, default: True
            Indicates whether the bins include the rightmost edge or not. If
            right == True (the default), then the bins [1,2,3,4] indicate
            (1,2], (2,3], (3,4].
        labels : array-like or bool, default: None
            Used as labels for the resulting bins. Must be of the same length as
            the resulting bins. If False, string bin labels are assigned by
            `pandas.cut`.
        precision : int, default: 3
            The precision at which to store and display the bins labels.
        include_lowest : bool, default: False
            Whether the first interval should be left-inclusive or not.
        squeeze : False
            This argument is deprecated.
        restore_coord_dims : bool, default: False
            If True, also restore the dimension order of multi-dimensional
            coordinates.
        duplicates : {"raise", "drop"}, default: "raise"
            If bin edges are not unique, raise ValueError or drop non-uniques.
        eagerly_compute_group: False, optional
            This argument is deprecated.

        Returns
        -------
        grouped : DatasetGroupBy
            A `DatasetGroupBy` object patterned after `pandas.GroupBy` that can be
            iterated over in the form of `(unique_value, grouped_array)` pairs.
            The name of the group has the added suffix `_bins` in order to
            distinguish it from the original variable.

        See Also
        --------
        :ref:`groupby`
            Users guide explanation of how to group and bin data.
        Dataset.groupby
        DataArray.groupby_bins
        core.groupby.DatasetGroupBy
        pandas.DataFrame.groupby

        References
        ----------
        .. [1] https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.cut.html
        """
        from xarray.core.groupby import (
            DatasetGroupBy,
            ResolvedGrouper,
            _validate_groupby_squeeze,
        )
        from xarray.groupers import BinGrouper

        _validate_groupby_squeeze(squeeze)
        grouper = BinGrouper(
            bins=bins,
            right=right,
            labels=labels,
            precision=precision,
            include_lowest=include_lowest,
        )
        rgrouper = ResolvedGrouper(
            grouper, group, self, eagerly_compute_group=eagerly_compute_group
        )

        return DatasetGroupBy(
            self,
            (rgrouper,),
            restore_coord_dims=restore_coord_dims,
        )

    def weighted(self, weights: DataArray) -> DatasetWeighted:
        """
        Weighted Dataset operations.

        Parameters
        ----------
        weights : DataArray
            An array of weights associated with the values in this Dataset.
            Each value in the data contributes to the reduction operation
            according to its associated weight.

        Notes
        -----
        ``weights`` must be a DataArray and cannot contain missing values.
        Missing values can be replaced by ``weights.fillna(0)``.

        Returns
        -------
        computation.weighted.DatasetWeighted

        See Also
        --------
        :func:`DataArray.weighted <DataArray.weighted>`

        :ref:`compute.weighted`
            User guide on weighted array reduction using :py:func:`~xarray.Dataset.weighted`

        :doc:`xarray-tutorial:fundamentals/03.4_weighted`
            Tutorial on Weighted Reduction using :py:func:`~xarray.Dataset.weighted`

        """
        from xarray.computation.weighted import DatasetWeighted

        return DatasetWeighted(self, weights)

    def rolling(
        self,
        dim: Mapping[Any, int] | None = None,
        min_periods: int | None = None,
        center: bool | Mapping[Any, bool] = False,
        **window_kwargs: int,
    ) -> DatasetRolling:
        """
        Rolling window object for Datasets.

        Parameters
        ----------
        dim : dict, optional
            Mapping from the dimension name to create the rolling iterator
            along (e.g. `time`) to its moving window size.
        min_periods : int or None, default: None
            Minimum number of observations in window required to have a value
            (otherwise result is NA). The default, None, is equivalent to
            setting min_periods equal to the size of the window.
        center : bool or Mapping to int, default: False
            Set the labels at the center of the window. The default, False,
            sets the labels at the right edge of the window.
        **window_kwargs : optional
            The keyword arguments form of ``dim``.
            One of dim or window_kwargs must be provided.

        Returns
        -------
        computation.rolling.DatasetRolling

        See Also
        --------
        Dataset.cumulative
        DataArray.rolling
        DataArray.rolling_exp
        """
        from xarray.computation.rolling import DatasetRolling

        dim = either_dict_or_kwargs(dim, window_kwargs, "rolling")
        return DatasetRolling(self, dim, min_periods=min_periods, center=center)

    def cumulative(
        self,
        dim: str | Iterable[Hashable],
        min_periods: int = 1,
    ) -> DatasetRolling:
        """
        Accumulating object for Datasets

        Parameters
        ----------
        dims : iterable of hashable
            The name(s) of the dimensions to create the cumulative window along
        min_periods : int, default: 1
            Minimum number of observations in window required to have a value
            (otherwise result is NA). The default is 1 (note this is different
            from ``Rolling``, whose default is the size of the window).

        Returns
        -------
        computation.rolling.DatasetRolling

        See Also
        --------
        DataArray.cumulative
        Dataset.rolling
        Dataset.rolling_exp
        """
        from xarray.computation.rolling import DatasetRolling

        if isinstance(dim, str):
            if dim not in self.dims:
                raise ValueError(
                    f"Dimension {dim} not found in data dimensions: {self.dims}"
                )
            dim = {dim: self.sizes[dim]}
        else:
            missing_dims = set(dim) - set(self.dims)
            if missing_dims:
                raise ValueError(
                    f"Dimensions {missing_dims} not found in data dimensions: {self.dims}"
                )
            dim = {d: self.sizes[d] for d in dim}

        return DatasetRolling(self, dim, min_periods=min_periods, center=False)

    def coarsen(
        self,
        dim: Mapping[Any, int] | None = None,
        boundary: CoarsenBoundaryOptions = "exact",
        side: SideOptions | Mapping[Any, SideOptions] = "left",
        coord_func: str | Callable | Mapping[Any, str | Callable] = "mean",
        **window_kwargs: int,
    ) -> DatasetCoarsen:
        """
        Coarsen object for Datasets.

        Parameters
        ----------
        dim : mapping of hashable to int, optional
            Mapping from the dimension name to the window size.
        boundary : {"exact", "trim", "pad"}, default: "exact"
            If 'exact', a ValueError will be raised if dimension size is not a
            multiple of the window size. If 'trim', the excess entries are
            dropped. If 'pad', NA will be padded.
        side : {"left", "right"} or mapping of str to {"left", "right"}, default: "left"
        coord_func : str or mapping of hashable to str, default: "mean"
            function (name) that is applied to the coordinates,
            or a mapping from coordinate name to function (name).

        Returns
        -------
        computation.rolling.DatasetCoarsen

        See Also
        --------
        :class:`computation.rolling.DatasetCoarsen`
        :func:`DataArray.coarsen <DataArray.coarsen>`

        :ref:`reshape.coarsen`
            User guide describing :py:func:`~xarray.Dataset.coarsen`

        :ref:`compute.coarsen`
            User guide on block arrgragation :py:func:`~xarray.Dataset.coarsen`

        :doc:`xarray-tutorial:fundamentals/03.3_windowed`
            Tutorial on windowed computation using :py:func:`~xarray.Dataset.coarsen`

        """
        from xarray.computation.rolling import DatasetCoarsen

        dim = either_dict_or_kwargs(dim, window_kwargs, "coarsen")
        return DatasetCoarsen(
            self,
            dim,
            boundary=boundary,
            side=side,
            coord_func=coord_func,
        )

    @_deprecate_positional_args("v2024.07.0")
    def resample(
        self,
        indexer: Mapping[Any, ResampleCompatible | Resampler] | None = None,
        *,
        skipna: bool | None = None,
        closed: SideOptions | None = None,
        label: SideOptions | None = None,
        offset: pd.Timedelta | datetime.timedelta | str | None = None,
        origin: str | DatetimeLike = "start_day",
        restore_coord_dims: bool | None = None,
        **indexer_kwargs: ResampleCompatible | Resampler,
    ) -> DatasetResample:
        """Returns a Resample object for performing resampling operations.

        Handles both downsampling and upsampling. The resampled
        dimension must be a datetime-like coordinate. If any intervals
        contain no values from the original object, they will be given
        the value ``NaN``.

        Parameters
        ----------
        indexer : Mapping of Hashable to str, datetime.timedelta, pd.Timedelta, pd.DateOffset, or Resampler, optional
            Mapping from the dimension name to resample frequency [1]_. The
            dimension must be datetime-like.
        skipna : bool, optional
            Whether to skip missing values when aggregating in downsampling.
        closed : {"left", "right"}, optional
            Side of each interval to treat as closed.
        label : {"left", "right"}, optional
            Side of each interval to use for labeling.
        origin : {'epoch', 'start', 'start_day', 'end', 'end_day'}, pd.Timestamp, datetime.datetime, np.datetime64, or cftime.datetime, default 'start_day'
            The datetime on which to adjust the grouping. The timezone of origin
            must match the timezone of the index.

            If a datetime is not used, these values are also supported:
            - 'epoch': `origin` is 1970-01-01
            - 'start': `origin` is the first value of the timeseries
            - 'start_day': `origin` is the first day at midnight of the timeseries
            - 'end': `origin` is the last value of the timeseries
            - 'end_day': `origin` is the ceiling midnight of the last day
        offset : pd.Timedelta, datetime.timedelta, or str, default is None
            An offset timedelta added to the origin.
        restore_coord_dims : bool, optional
            If True, also restore the dimension order of multi-dimensional
            coordinates.
        **indexer_kwargs : str, datetime.timedelta, pd.Timedelta, pd.DateOffset, or Resampler
            The keyword arguments form of ``indexer``.
            One of indexer or indexer_kwargs must be provided.

        Returns
        -------
        resampled : core.resample.DataArrayResample
            This object resampled.

        See Also
        --------
        DataArray.resample
        pandas.Series.resample
        pandas.DataFrame.resample
        Dataset.groupby
        DataArray.groupby

        References
        ----------
        .. [1] https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases
        """
        from xarray.core.resample import DatasetResample

        return self._resample(
            resample_cls=DatasetResample,
            indexer=indexer,
            skipna=skipna,
            closed=closed,
            label=label,
            offset=offset,
            origin=origin,
            restore_coord_dims=restore_coord_dims,
            **indexer_kwargs,
        )

    def drop_attrs(self, *, deep: bool = True) -> Self:
        """
        Removes all attributes from the Dataset and its variables.

        Parameters
        ----------
        deep : bool, default True
            Removes attributes from all variables.

        Returns
        -------
        Dataset
        """
        # Remove attributes from the dataset
        self = self._replace(attrs={})

        if not deep:
            return self

        # Remove attributes from each variable in the dataset
        for var in self.variables:
            # variables don't have a `._replace` method, so we copy and then remove
            # attrs. If we added a `._replace` method, we could use that instead.
            if var not in self.indexes:
                self[var] = self[var].copy()
                self[var].attrs = {}

        new_idx_variables = {}
        # Not sure this is the most elegant way of doing this, but it works.
        # (Should we have a more general "map over all variables, including
        # indexes" approach?)
        for idx, idx_vars in self.xindexes.group_by_index():
            # copy each coordinate variable of an index and drop their attrs
            temp_idx_variables = {k: v.copy() for k, v in idx_vars.items()}
            for v in temp_idx_variables.values():
                v.attrs = {}
            # re-wrap the index object in new coordinate variables
            new_idx_variables.update(idx.create_variables(temp_idx_variables))
        self = self.assign(new_idx_variables)

        return self
