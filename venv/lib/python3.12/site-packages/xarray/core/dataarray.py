from __future__ import annotations

import copy
import datetime
import warnings
from collections.abc import (
    Callable,
    Hashable,
    Iterable,
    Mapping,
    MutableMapping,
    Sequence,
)
from functools import partial
from os import PathLike
from types import EllipsisType
from typing import TYPE_CHECKING, Any, Generic, Literal, NoReturn, TypeVar, overload

import numpy as np
import pandas as pd

from xarray.coding.calendar_ops import convert_calendar, interp_calendar
from xarray.coding.cftimeindex import CFTimeIndex
from xarray.computation import computation, ops
from xarray.computation.arithmetic import DataArrayArithmetic
from xarray.core import dtypes, indexing, utils
from xarray.core._aggregations import DataArrayAggregations
from xarray.core.accessor_dt import CombinedDatetimelikeAccessor
from xarray.core.accessor_str import StringAccessor
from xarray.core.common import AbstractArray, DataWithCoords, get_chunksizes
from xarray.core.coordinates import (
    Coordinates,
    DataArrayCoordinates,
    assert_coordinate_consistent,
    create_coords_with_default_indexes,
    validate_dataarray_coords,
)
from xarray.core.dataset import Dataset
from xarray.core.extension_array import PandasExtensionArray
from xarray.core.formatting import format_item
from xarray.core.indexes import (
    Index,
    Indexes,
    PandasMultiIndex,
    filter_indexes_from_coords,
    isel_indexes,
)
from xarray.core.indexing import is_fancy_indexer, map_index_queries
from xarray.core.options import OPTIONS, _get_keep_attrs
from xarray.core.types import (
    Bins,
    DaCompatible,
    NetcdfWriteModes,
    T_Chunks,
    T_DataArray,
    T_DataArrayOrSet,
    ZarrWriteModes,
)
from xarray.core.utils import (
    Default,
    FilteredMapping,
    ReprObject,
    _default,
    either_dict_or_kwargs,
    hashable,
    infix_dims,
    result_name,
)
from xarray.core.variable import (
    IndexVariable,
    Variable,
    as_compatible_data,
    as_variable,
)
from xarray.plot.accessor import DataArrayPlotAccessor
from xarray.plot.utils import _get_units_from_attrs
from xarray.structure import alignment
from xarray.structure.alignment import (
    _broadcast_helper,
    _get_broadcast_dims_map_common_coords,
    align,
)
from xarray.structure.chunks import unify_chunks
from xarray.structure.merge import PANDAS_TYPES, MergeError
from xarray.util.deprecation_helpers import _deprecate_positional_args, deprecate_dims

if TYPE_CHECKING:
    from dask.dataframe import DataFrame as DaskDataFrame
    from dask.delayed import Delayed
    from iris.cube import Cube as iris_Cube
    from numpy.typing import ArrayLike

    from xarray.backends import ZarrStore
    from xarray.backends.api import T_NetcdfEngine, T_NetcdfTypes
    from xarray.computation.rolling import DataArrayCoarsen, DataArrayRolling
    from xarray.computation.weighted import DataArrayWeighted
    from xarray.core.groupby import DataArrayGroupBy
    from xarray.core.resample import DataArrayResample
    from xarray.core.types import (
        CoarsenBoundaryOptions,
        DatetimeLike,
        DatetimeUnitOptions,
        Dims,
        ErrorOptions,
        ErrorOptionsWithWarn,
        GroupIndices,
        GroupInput,
        InterpOptions,
        PadModeOptions,
        PadReflectOptions,
        QuantileMethods,
        QueryEngineOptions,
        QueryParserOptions,
        ReindexMethodOptions,
        ResampleCompatible,
        Self,
        SideOptions,
        T_ChunkDimFreq,
        T_ChunksFreq,
        T_Xarray,
    )
    from xarray.groupers import Grouper, Resampler
    from xarray.namedarray.parallelcompat import ChunkManagerEntrypoint

    T_XarrayOther = TypeVar("T_XarrayOther", bound="DataArray" | Dataset)


def _infer_coords_and_dims(
    shape: tuple[int, ...],
    coords: (
        Sequence[Sequence | pd.Index | DataArray | Variable | np.ndarray]
        | Mapping
        | None
    ),
    dims: str | Iterable[Hashable] | None,
) -> tuple[Mapping[Hashable, Any], tuple[Hashable, ...]]:
    """All the logic for creating a new DataArray"""

    if (
        coords is not None
        and not utils.is_dict_like(coords)
        and len(coords) != len(shape)
    ):
        raise ValueError(
            f"coords is not dict-like, but it has {len(coords)} items, "
            f"which does not match the {len(shape)} dimensions of the "
            "data"
        )

    if isinstance(dims, str):
        dims = (dims,)
    elif dims is None:
        dims = [f"dim_{n}" for n in range(len(shape))]
        if coords is not None and len(coords) == len(shape):
            # try to infer dimensions from coords
            if utils.is_dict_like(coords):
                dims = list(coords.keys())
            else:
                for n, (dim, coord) in enumerate(zip(dims, coords, strict=True)):
                    coord = as_variable(
                        coord, name=dim, auto_convert=False
                    ).to_index_variable()
                    dims[n] = coord.name
    dims_tuple = tuple(dims)
    if len(dims_tuple) != len(shape):
        raise ValueError(
            "different number of dimensions on data "
            f"and dims: {len(shape)} vs {len(dims_tuple)}"
        )
    for d in dims_tuple:
        if not hashable(d):
            raise TypeError(f"Dimension {d} is not hashable")

    new_coords: Mapping[Hashable, Any]

    if isinstance(coords, Coordinates):
        new_coords = coords
    else:
        new_coords = {}
        if utils.is_dict_like(coords):
            for k, v in coords.items():
                new_coords[k] = as_variable(v, name=k, auto_convert=False)
                if new_coords[k].dims == (k,):
                    new_coords[k] = new_coords[k].to_index_variable()
        elif coords is not None:
            for dim, coord in zip(dims_tuple, coords, strict=True):
                var = as_variable(coord, name=dim, auto_convert=False)
                var.dims = (dim,)
                new_coords[dim] = var.to_index_variable()

    validate_dataarray_coords(shape, new_coords, dims_tuple)

    return new_coords, dims_tuple


def _check_data_shape(
    data: Any,
    coords: (
        Sequence[Sequence | pd.Index | DataArray | Variable | np.ndarray]
        | Mapping
        | None
    ),
    dims: str | Iterable[Hashable] | None,
) -> Any:
    if data is dtypes.NA:
        data = np.nan
    if coords is not None and utils.is_scalar(data, include_0d=False):
        if utils.is_dict_like(coords):
            if dims is None:
                return data
            else:
                data_shape = tuple(
                    (
                        as_variable(coords[k], k, auto_convert=False).size
                        if k in coords.keys()
                        else 1
                    )
                    for k in dims
                )
        else:
            data_shape = tuple(
                as_variable(coord, "foo", auto_convert=False).size for coord in coords
            )
        data = np.full(data_shape, data)
    return data


class _LocIndexer(Generic[T_DataArray]):
    __slots__ = ("data_array",)

    def __init__(self, data_array: T_DataArray):
        self.data_array = data_array

    def __getitem__(self, key) -> T_DataArray:
        if not utils.is_dict_like(key):
            # expand the indexer so we can handle Ellipsis
            labels = indexing.expanded_indexer(key, self.data_array.ndim)
            key = dict(zip(self.data_array.dims, labels, strict=True))
        return self.data_array.sel(key)

    def __setitem__(self, key, value) -> None:
        if not utils.is_dict_like(key):
            # expand the indexer so we can handle Ellipsis
            labels = indexing.expanded_indexer(key, self.data_array.ndim)
            key = dict(zip(self.data_array.dims, labels, strict=True))

        dim_indexers = map_index_queries(self.data_array, key).dim_indexers
        self.data_array[dim_indexers] = value


# Used as the key corresponding to a DataArray's variable when converting
# arbitrary DataArray objects to datasets
_THIS_ARRAY = ReprObject("<this-array>")


class DataArray(
    AbstractArray,
    DataWithCoords,
    DataArrayArithmetic,
    DataArrayAggregations,
):
    """N-dimensional array with labeled coordinates and dimensions.

    DataArray provides a wrapper around numpy ndarrays that uses
    labeled dimensions and coordinates to support metadata aware
    operations. The API is similar to that for the pandas Series or
    DataFrame, but DataArray objects can have any number of dimensions,
    and their contents have fixed data types.

    Additional features over raw numpy arrays:

    - Apply operations over dimensions by name: ``x.sum('time')``.
    - Select or assign values by integer location (like numpy):
      ``x[:10]`` or by label (like pandas): ``x.loc['2014-01-01']`` or
      ``x.sel(time='2014-01-01')``.
    - Mathematical operations (e.g., ``x - y``) vectorize across
      multiple dimensions (known in numpy as "broadcasting") based on
      dimension names, regardless of their original order.
    - Keep track of arbitrary metadata in the form of a Python
      dictionary: ``x.attrs``
    - Convert to a pandas Series: ``x.to_series()``.

    Getting items from or doing mathematical operations with a
    DataArray always returns another DataArray.

    Parameters
    ----------
    data : array_like
        Values for this array. Must be an ``numpy.ndarray``, ndarray
        like, or castable to an ``ndarray``. If a self-described xarray
        or pandas object, attempts are made to use this array's
        metadata to fill in other unspecified arguments. A view of the
        array's data is used instead of a copy if possible.
    coords : sequence or dict of array_like or :py:class:`~xarray.Coordinates`, optional
        Coordinates (tick labels) to use for indexing along each
        dimension. The following notations are accepted:

        - mapping {dimension name: array-like}
        - sequence of tuples that are valid arguments for
          ``xarray.Variable()``
          - (dims, data)
          - (dims, data, attrs)
          - (dims, data, attrs, encoding)

        Additionally, it is possible to define a coord whose name
        does not match the dimension name, or a coord based on multiple
        dimensions, with one of the following notations:

        - mapping {coord name: DataArray}
        - mapping {coord name: Variable}
        - mapping {coord name: (dimension name, array-like)}
        - mapping {coord name: (tuple of dimension names, array-like)}

        Alternatively, a :py:class:`~xarray.Coordinates` object may be used in
        order to explicitly pass indexes (e.g., a multi-index or any custom
        Xarray index) or to bypass the creation of a default index for any
        :term:`Dimension coordinate` included in that object.
    dims : Hashable or sequence of Hashable, optional
        Name(s) of the data dimension(s). Must be either a Hashable
        (only for 1D data) or a sequence of Hashables with length equal
        to the number of dimensions. If this argument is omitted,
        dimension names are taken from ``coords`` (if possible) and
        otherwise default to ``['dim_0', ... 'dim_n']``.
    name : str or None, optional
        Name of this array.
    attrs : dict_like or None, optional
        Attributes to assign to the new instance. By default, an empty
        attribute dictionary is initialized.
        (see FAQ, :ref:`approach to metadata`)
    indexes : :py:class:`~xarray.Indexes` or dict-like, optional
        For internal use only. For passing indexes objects to the
        new DataArray, use the ``coords`` argument instead with a
        :py:class:`~xarray.Coordinate` object (both coordinate variables
        and indexes will be extracted from the latter).

    Examples
    --------
    Create data:

    >>> np.random.seed(0)
    >>> temperature = 15 + 8 * np.random.randn(2, 2, 3)
    >>> lon = [[-99.83, -99.32], [-99.79, -99.23]]
    >>> lat = [[42.25, 42.21], [42.63, 42.59]]
    >>> time = pd.date_range("2014-09-06", periods=3)
    >>> reference_time = pd.Timestamp("2014-09-05")

    Initialize a dataarray with multiple dimensions:

    >>> da = xr.DataArray(
    ...     data=temperature,
    ...     dims=["x", "y", "time"],
    ...     coords=dict(
    ...         lon=(["x", "y"], lon),
    ...         lat=(["x", "y"], lat),
    ...         time=time,
    ...         reference_time=reference_time,
    ...     ),
    ...     attrs=dict(
    ...         description="Ambient temperature.",
    ...         units="degC",
    ...     ),
    ... )
    >>> da
    <xarray.DataArray (x: 2, y: 2, time: 3)> Size: 96B
    array([[[29.11241877, 18.20125767, 22.82990387],
            [32.92714559, 29.94046392,  7.18177696]],
    <BLANKLINE>
           [[22.60070734, 13.78914233, 14.17424919],
            [18.28478802, 16.15234857, 26.63418806]]])
    Coordinates:
        lon             (x, y) float64 32B -99.83 -99.32 -99.79 -99.23
        lat             (x, y) float64 32B 42.25 42.21 42.63 42.59
      * time            (time) datetime64[ns] 24B 2014-09-06 2014-09-07 2014-09-08
        reference_time  datetime64[ns] 8B 2014-09-05
    Dimensions without coordinates: x, y
    Attributes:
        description:  Ambient temperature.
        units:        degC

    Find out where the coldest temperature was:

    >>> da.isel(da.argmin(...))
    <xarray.DataArray ()> Size: 8B
    array(7.18177696)
    Coordinates:
        lon             float64 8B -99.32
        lat             float64 8B 42.21
        time            datetime64[ns] 8B 2014-09-08
        reference_time  datetime64[ns] 8B 2014-09-05
    Attributes:
        description:  Ambient temperature.
        units:        degC
    """

    _cache: dict[str, Any]
    _coords: dict[Any, Variable]
    _close: Callable[[], None] | None
    _indexes: dict[Hashable, Index]
    _name: Hashable | None
    _variable: Variable

    __slots__ = (
        "__weakref__",
        "_cache",
        "_close",
        "_coords",
        "_indexes",
        "_name",
        "_variable",
    )

    dt = utils.UncachedAccessor(CombinedDatetimelikeAccessor["DataArray"])

    def __init__(
        self,
        data: Any = dtypes.NA,
        coords: (
            Sequence[Sequence | pd.Index | DataArray | Variable | np.ndarray]
            | Mapping
            | None
        ) = None,
        dims: str | Iterable[Hashable] | None = None,
        name: Hashable | None = None,
        attrs: Mapping | None = None,
        # internal parameters
        indexes: Mapping[Hashable, Index] | None = None,
        fastpath: bool = False,
    ) -> None:
        if fastpath:
            variable = data
            assert dims is None
            assert attrs is None
            assert indexes is not None
        else:
            if indexes is not None:
                raise ValueError(
                    "Explicitly passing indexes via the `indexes` argument is not supported "
                    "when `fastpath=False`. Use the `coords` argument instead."
                )

            # try to fill in arguments from data if they weren't supplied
            if coords is None:
                if isinstance(data, DataArray):
                    coords = data.coords
                elif isinstance(data, pd.Series):
                    coords = [data.index]
                elif isinstance(data, pd.DataFrame):
                    coords = [data.index, data.columns]
                elif isinstance(data, pd.Index | IndexVariable):
                    coords = [data]

            if dims is None:
                dims = getattr(data, "dims", getattr(coords, "dims", None))
            if name is None:
                name = getattr(data, "name", None)
            if attrs is None and not isinstance(data, PANDAS_TYPES):
                attrs = getattr(data, "attrs", None)

            data = _check_data_shape(data, coords, dims)
            data = as_compatible_data(data)
            coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
            variable = Variable(dims, data, attrs, fastpath=True)

            if not isinstance(coords, Coordinates):
                coords = create_coords_with_default_indexes(coords)
            indexes = dict(coords.xindexes)
            coords = {k: v.copy() for k, v in coords.variables.items()}

        # These fully describe a DataArray
        self._variable = variable
        assert isinstance(coords, dict)
        self._coords = coords
        self._name = name
        self._indexes = dict(indexes)

        self._close = None

    @classmethod
    def _construct_direct(
        cls,
        variable: Variable,
        coords: dict[Any, Variable],
        name: Hashable,
        indexes: dict[Hashable, Index],
    ) -> Self:
        """Shortcut around __init__ for internal use when we want to skip
        costly validation
        """
        obj = object.__new__(cls)
        obj._variable = variable
        obj._coords = coords
        obj._name = name
        obj._indexes = indexes
        obj._close = None
        return obj

    def _replace(
        self,
        variable: Variable | None = None,
        coords=None,
        name: Hashable | Default | None = _default,
        attrs=_default,
        indexes=None,
    ) -> Self:
        if variable is None:
            variable = self.variable
        if coords is None:
            coords = self._coords
        if indexes is None:
            indexes = self._indexes
        if name is _default:
            name = self.name
        if attrs is _default:
            attrs = copy.copy(self.attrs)
        else:
            variable = variable.copy()
            variable.attrs = attrs
        return type(self)(variable, coords, name=name, indexes=indexes, fastpath=True)

    def _replace_maybe_drop_dims(
        self,
        variable: Variable,
        name: Hashable | Default | None = _default,
    ) -> Self:
        if self.sizes == variable.sizes:
            coords = self._coords.copy()
            indexes = self._indexes
        elif set(self.dims) == set(variable.dims):
            # Shape has changed (e.g. from reduce(..., keepdims=True)
            new_sizes = dict(zip(self.dims, variable.shape, strict=True))
            coords = {
                k: v
                for k, v in self._coords.items()
                if v.shape == tuple(new_sizes[d] for d in v.dims)
            }
            indexes = filter_indexes_from_coords(self._indexes, set(coords))
        else:
            allowed_dims = set(variable.dims)
            coords = {
                k: v for k, v in self._coords.items() if set(v.dims) <= allowed_dims
            }
            indexes = filter_indexes_from_coords(self._indexes, set(coords))
        return self._replace(variable, coords, name, indexes=indexes)

    def _overwrite_indexes(
        self,
        indexes: Mapping[Any, Index],
        variables: Mapping[Any, Variable] | None = None,
        drop_coords: list[Hashable] | None = None,
        rename_dims: Mapping[Any, Any] | None = None,
    ) -> Self:
        """Maybe replace indexes and their corresponding coordinates."""
        if not indexes:
            return self

        if variables is None:
            variables = {}
        if drop_coords is None:
            drop_coords = []

        new_variable = self.variable.copy()
        new_coords = self._coords.copy()
        new_indexes = dict(self._indexes)

        for name in indexes:
            new_coords[name] = variables[name]
            new_indexes[name] = indexes[name]

        for name in drop_coords:
            new_coords.pop(name)
            new_indexes.pop(name)

        if rename_dims:
            new_variable.dims = tuple(rename_dims.get(d, d) for d in new_variable.dims)

        return self._replace(
            variable=new_variable, coords=new_coords, indexes=new_indexes
        )

    def _to_temp_dataset(self) -> Dataset:
        return self._to_dataset_whole(name=_THIS_ARRAY, shallow_copy=False)

    def _from_temp_dataset(
        self, dataset: Dataset, name: Hashable | Default | None = _default
    ) -> Self:
        variable = dataset._variables.pop(_THIS_ARRAY)
        coords = dataset._variables
        indexes = dataset._indexes
        return self._replace(variable, coords, name, indexes=indexes)

    def _to_dataset_split(self, dim: Hashable) -> Dataset:
        """splits dataarray along dimension 'dim'"""

        def subset(dim, label):
            array = self.loc[{dim: label}]
            array.attrs = {}
            return as_variable(array)

        variables_from_split = {
            label: subset(dim, label) for label in self.get_index(dim)
        }
        coord_names = set(self._coords) - {dim}

        ambiguous_vars = set(variables_from_split) & coord_names
        if ambiguous_vars:
            rename_msg_fmt = ", ".join([f"{v}=..." for v in sorted(ambiguous_vars)])
            raise ValueError(
                f"Splitting along the dimension {dim!r} would produce the variables "
                f"{tuple(sorted(ambiguous_vars))} which are also existing coordinate "
                f"variables. Use DataArray.rename({rename_msg_fmt}) or "
                f"DataArray.assign_coords({dim}=...) to resolve this ambiguity."
            )

        variables = variables_from_split | {
            k: v for k, v in self._coords.items() if k != dim
        }
        indexes = filter_indexes_from_coords(self._indexes, coord_names)
        dataset = Dataset._construct_direct(
            variables, coord_names, indexes=indexes, attrs=self.attrs
        )
        return dataset

    def _to_dataset_whole(
        self, name: Hashable = None, shallow_copy: bool = True
    ) -> Dataset:
        if name is None:
            name = self.name
        if name is None:
            raise ValueError(
                "unable to convert unnamed DataArray to a "
                "Dataset without providing an explicit name"
            )
        if name in self.coords:
            raise ValueError(
                "cannot create a Dataset from a DataArray with "
                "the same name as one of its coordinates"
            )
        # use private APIs for speed: this is called by _to_temp_dataset(),
        # which is used in the guts of a lot of operations (e.g., reindex)
        variables = self._coords.copy()
        variables[name] = self.variable
        if shallow_copy:
            for k in variables:
                variables[k] = variables[k].copy(deep=False)
        indexes = self._indexes

        coord_names = set(self._coords)
        return Dataset._construct_direct(variables, coord_names, indexes=indexes)

    def to_dataset(
        self,
        dim: Hashable = None,
        *,
        name: Hashable = None,
        promote_attrs: bool = False,
    ) -> Dataset:
        """Convert a DataArray to a Dataset.

        Parameters
        ----------
        dim : Hashable, optional
            Name of the dimension on this array along which to split this array
            into separate variables. If not provided, this array is converted
            into a Dataset of one variable.
        name : Hashable, optional
            Name to substitute for this array's name. Only valid if ``dim`` is
            not provided.
        promote_attrs : bool, default: False
            Set to True to shallow copy attrs of DataArray to returned Dataset.

        Returns
        -------
        dataset : Dataset
        """
        if dim is not None and dim not in self.dims:
            raise TypeError(
                f"{dim} is not a dim. If supplying a ``name``, pass as a kwarg."
            )

        if dim is not None:
            if name is not None:
                raise TypeError("cannot supply both dim and name arguments")
            result = self._to_dataset_split(dim)
        else:
            result = self._to_dataset_whole(name)

        if promote_attrs:
            result.attrs = dict(self.attrs)

        return result

    @property
    def name(self) -> Hashable | None:
        """The name of this array."""
        return self._name

    @name.setter
    def name(self, value: Hashable | None) -> None:
        self._name = value

    @property
    def variable(self) -> Variable:
        """Low level interface to the Variable object for this DataArray."""
        return self._variable

    @property
    def dtype(self) -> np.dtype:
        """
        Data-type of the array’s elements.

        See Also
        --------
        ndarray.dtype
        numpy.dtype
        """
        return self.variable.dtype

    @property
    def shape(self) -> tuple[int, ...]:
        """
        Tuple of array dimensions.

        See Also
        --------
        numpy.ndarray.shape
        """
        return self.variable.shape

    @property
    def size(self) -> int:
        """
        Number of elements in the array.

        Equal to ``np.prod(a.shape)``, i.e., the product of the array’s dimensions.

        See Also
        --------
        numpy.ndarray.size
        """
        return self.variable.size

    @property
    def nbytes(self) -> int:
        """
        Total bytes consumed by the elements of this DataArray's data.

        If the underlying data array does not include ``nbytes``, estimates
        the bytes consumed based on the ``size`` and ``dtype``.
        """
        return self.variable.nbytes

    @property
    def ndim(self) -> int:
        """
        Number of array dimensions.

        See Also
        --------
        numpy.ndarray.ndim
        """
        return self.variable.ndim

    def __len__(self) -> int:
        return len(self.variable)

    @property
    def data(self) -> Any:
        """
        The DataArray's data as an array. The underlying array type
        (e.g. dask, sparse, pint) is preserved.

        See Also
        --------
        DataArray.to_numpy
        DataArray.as_numpy
        DataArray.values
        """
        return self.variable.data

    @data.setter
    def data(self, value: Any) -> None:
        self.variable.data = value

    @property
    def values(self) -> np.ndarray:
        """
        The array's data converted to numpy.ndarray.

        This will attempt to convert the array naively using np.array(),
        which will raise an error if the array type does not support
        coercion like this (e.g. cupy).

        Note that this array is not copied; operations on it follow
        numpy's rules of what generates a view vs. a copy, and changes
        to this array may be reflected in the DataArray as well.
        """
        return self.variable.values

    @values.setter
    def values(self, value: Any) -> None:
        self.variable.values = value

    def to_numpy(self) -> np.ndarray:
        """
        Coerces wrapped data to numpy and returns a numpy.ndarray.

        See Also
        --------
        DataArray.as_numpy : Same but returns the surrounding DataArray instead.
        Dataset.as_numpy
        DataArray.values
        DataArray.data
        """
        return self.variable.to_numpy()

    def as_numpy(self) -> Self:
        """
        Coerces wrapped data and coordinates into numpy arrays, returning a DataArray.

        See Also
        --------
        DataArray.to_numpy : Same but returns only the data as a numpy.ndarray object.
        Dataset.as_numpy : Converts all variables in a Dataset.
        DataArray.values
        DataArray.data
        """
        coords = {k: v.as_numpy() for k, v in self._coords.items()}
        return self._replace(self.variable.as_numpy(), coords, indexes=self._indexes)

    @property
    def _in_memory(self) -> bool:
        return self.variable._in_memory

    def _to_index(self) -> pd.Index:
        return self.variable._to_index()

    def to_index(self) -> pd.Index:
        """Convert this variable to a pandas.Index. Only possible for 1D
        arrays.
        """
        return self.variable.to_index()

    @property
    def dims(self) -> tuple[Hashable, ...]:
        """Tuple of dimension names associated with this array.

        Note that the type of this property is inconsistent with
        `Dataset.dims`.  See `Dataset.sizes` and `DataArray.sizes` for
        consistently named properties.

        See Also
        --------
        DataArray.sizes
        Dataset.dims
        """
        return self.variable.dims

    @dims.setter
    def dims(self, value: Any) -> NoReturn:
        raise AttributeError(
            "you cannot assign dims on a DataArray. Use "
            ".rename() or .swap_dims() instead."
        )

    def _item_key_to_dict(self, key: Any) -> Mapping[Hashable, Any]:
        if utils.is_dict_like(key):
            return key
        key = indexing.expanded_indexer(key, self.ndim)
        return dict(zip(self.dims, key, strict=True))

    def _getitem_coord(self, key: Any) -> Self:
        from xarray.core.dataset_utils import _get_virtual_variable

        try:
            var = self._coords[key]
        except KeyError:
            dim_sizes = dict(zip(self.dims, self.shape, strict=True))
            _, key, var = _get_virtual_variable(self._coords, key, dim_sizes)

        return self._replace_maybe_drop_dims(var, name=key)

    def __getitem__(self, key: Any) -> Self:
        if isinstance(key, str):
            return self._getitem_coord(key)
        else:
            # xarray-style array indexing
            return self.isel(indexers=self._item_key_to_dict(key))

    def __setitem__(self, key: Any, value: Any) -> None:
        if isinstance(key, str):
            self.coords[key] = value
        else:
            # Coordinates in key, value and self[key] should be consistent.
            # TODO Coordinate consistency in key is checked here, but it
            # causes unnecessary indexing. It should be optimized.
            obj = self[key]
            if isinstance(value, DataArray):
                assert_coordinate_consistent(value, obj.coords.variables)
                value = value.variable
            # DataArray key -> Variable key
            key = {
                k: v.variable if isinstance(v, DataArray) else v
                for k, v in self._item_key_to_dict(key).items()
            }
            self.variable[key] = value

    def __delitem__(self, key: Any) -> None:
        del self.coords[key]

    @property
    def _attr_sources(self) -> Iterable[Mapping[Hashable, Any]]:
        """Places to look-up items for attribute-style access"""
        yield from self._item_sources
        yield self.attrs

    @property
    def _item_sources(self) -> Iterable[Mapping[Hashable, Any]]:
        """Places to look-up items for key-completion"""
        yield FilteredMapping(keys=self._coords, mapping=self.coords)

        # virtual coordinates
        yield FilteredMapping(keys=self.dims, mapping=self.coords)

    def __contains__(self, key: Any) -> bool:
        return key in self.data

    @property
    def loc(self) -> _LocIndexer:
        """Attribute for location based indexing like pandas."""
        return _LocIndexer(self)

    @property
    def attrs(self) -> dict[Any, Any]:
        """Dictionary storing arbitrary metadata with this array."""
        return self.variable.attrs

    @attrs.setter
    def attrs(self, value: Mapping[Any, Any]) -> None:
        self.variable.attrs = dict(value)

    @property
    def encoding(self) -> dict[Any, Any]:
        """Dictionary of format-specific settings for how this array should be
        serialized."""
        return self.variable.encoding

    @encoding.setter
    def encoding(self, value: Mapping[Any, Any]) -> None:
        self.variable.encoding = dict(value)

    def reset_encoding(self) -> Self:
        warnings.warn(
            "reset_encoding is deprecated since 2023.11, use `drop_encoding` instead",
            stacklevel=2,
        )
        return self.drop_encoding()

    def drop_encoding(self) -> Self:
        """Return a new DataArray without encoding on the array or any attached
        coords."""
        ds = self._to_temp_dataset().drop_encoding()
        return self._from_temp_dataset(ds)

    @property
    def indexes(self) -> Indexes:
        """Mapping of pandas.Index objects used for label based indexing.

        Raises an error if this Dataset has indexes that cannot be coerced
        to pandas.Index objects.

        See Also
        --------
        DataArray.xindexes

        """
        return self.xindexes.to_pandas_indexes()

    @property
    def xindexes(self) -> Indexes[Index]:
        """Mapping of :py:class:`~xarray.indexes.Index` objects
        used for label based indexing.
        """
        return Indexes(self._indexes, {k: self._coords[k] for k in self._indexes})

    @property
    def coords(self) -> DataArrayCoordinates:
        """Mapping of :py:class:`~xarray.DataArray` objects corresponding to
        coordinate variables.

        See Also
        --------
        Coordinates
        """
        return DataArrayCoordinates(self)

    @overload
    def reset_coords(
        self,
        names: Dims = None,
        *,
        drop: Literal[False] = False,
    ) -> Dataset: ...

    @overload
    def reset_coords(
        self,
        names: Dims = None,
        *,
        drop: Literal[True],
    ) -> Self: ...

    def reset_coords(
        self,
        names: Dims = None,
        *,
        drop: bool = False,
    ) -> Self | Dataset:
        """Given names of coordinates, reset them to become variables.

        Parameters
        ----------
        names : str, Iterable of Hashable or None, optional
            Name(s) of non-index coordinates in this dataset to reset into
            variables. By default, all non-index coordinates are reset.
        drop : bool, default: False
            If True, remove coordinates instead of converting them into
            variables.

        Returns
        -------
        Dataset, or DataArray if ``drop == True``

        Examples
        --------
        >>> temperature = np.arange(25).reshape(5, 5)
        >>> pressure = np.arange(50, 75).reshape(5, 5)
        >>> da = xr.DataArray(
        ...     data=temperature,
        ...     dims=["x", "y"],
        ...     coords=dict(
        ...         lon=("x", np.arange(10, 15)),
        ...         lat=("y", np.arange(20, 25)),
        ...         Pressure=(["x", "y"], pressure),
        ...     ),
        ...     name="Temperature",
        ... )
        >>> da
        <xarray.DataArray 'Temperature' (x: 5, y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Coordinates:
            lon       (x) int64 40B 10 11 12 13 14
            lat       (y) int64 40B 20 21 22 23 24
            Pressure  (x, y) int64 200B 50 51 52 53 54 55 56 57 ... 68 69 70 71 72 73 74
        Dimensions without coordinates: x, y

        Return Dataset with target coordinate as a data variable rather than a coordinate variable:

        >>> da.reset_coords(names="Pressure")
        <xarray.Dataset> Size: 480B
        Dimensions:      (x: 5, y: 5)
        Coordinates:
            lon          (x) int64 40B 10 11 12 13 14
            lat          (y) int64 40B 20 21 22 23 24
        Dimensions without coordinates: x, y
        Data variables:
            Pressure     (x, y) int64 200B 50 51 52 53 54 55 56 ... 68 69 70 71 72 73 74
            Temperature  (x, y) int64 200B 0 1 2 3 4 5 6 7 8 ... 17 18 19 20 21 22 23 24

        Return DataArray without targeted coordinate:

        >>> da.reset_coords(names="Pressure", drop=True)
        <xarray.DataArray 'Temperature' (x: 5, y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Coordinates:
            lon      (x) int64 40B 10 11 12 13 14
            lat      (y) int64 40B 20 21 22 23 24
        Dimensions without coordinates: x, y
        """
        if names is None:
            names = set(self.coords) - set(self._indexes)
        dataset = self.coords.to_dataset().reset_coords(names, drop)
        if drop:
            return self._replace(coords=dataset._variables)
        if self.name is None:
            raise ValueError(
                "cannot reset_coords with drop=False on an unnamed DataArray"
            )
        dataset[self.name] = self.variable
        return dataset

    def __dask_tokenize__(self) -> object:
        from dask.base import normalize_token

        return normalize_token((type(self), self._variable, self._coords, self._name))

    def __dask_graph__(self):
        return self._to_temp_dataset().__dask_graph__()

    def __dask_keys__(self):
        return self._to_temp_dataset().__dask_keys__()

    def __dask_layers__(self):
        return self._to_temp_dataset().__dask_layers__()

    @property
    def __dask_optimize__(self):
        return self._to_temp_dataset().__dask_optimize__

    @property
    def __dask_scheduler__(self):
        return self._to_temp_dataset().__dask_scheduler__

    def __dask_postcompute__(self):
        func, args = self._to_temp_dataset().__dask_postcompute__()
        return self._dask_finalize, (self.name, func) + args

    def __dask_postpersist__(self):
        func, args = self._to_temp_dataset().__dask_postpersist__()
        return self._dask_finalize, (self.name, func) + args

    @classmethod
    def _dask_finalize(cls, results, name, func, *args, **kwargs) -> Self:
        ds = func(results, *args, **kwargs)
        variable = ds._variables.pop(_THIS_ARRAY)
        coords = ds._variables
        indexes = ds._indexes
        return cls(variable, coords, name=name, indexes=indexes, fastpath=True)

    def load(self, **kwargs) -> Self:
        """Trigger loading data into memory and return this dataarray.

        Data will be computed and/or loaded from disk or a remote source.

        Unlike ``.compute``, the original dataarray is modified and returned.

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
        object : DataArray
            Same object but with lazy data and coordinates as in-memory arrays.

        See Also
        --------
        dask.compute
        DataArray.load_async
        DataArray.compute
        Dataset.load
        Variable.load
        """
        ds = self._to_temp_dataset().load(**kwargs)
        new = self._from_temp_dataset(ds)
        self._variable = new._variable
        self._coords = new._coords
        return self

    async def load_async(self, **kwargs) -> Self:
        """Trigger and await asynchronous loading of data into memory and return this dataarray.

        Data will be computed and/or loaded from disk or a remote source.

        Unlike ``.compute``, the original dataarray is modified and returned.

        Only works when opening data lazily from IO storage backends which support lazy asynchronous loading.
        Otherwise will raise a NotImplementedError.

        Note users are expected to limit concurrency themselves - xarray does not internally limit concurrency in any way.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.compute``.

        Returns
        -------
        object : Dataarray
            Same object but with lazy data and coordinates as in-memory arrays.

        See Also
        --------
        dask.compute
        DataArray.compute
        DataArray.load
        Dataset.load_async
        Variable.load_async
        """
        temp_ds = self._to_temp_dataset()
        ds = await temp_ds.load_async(**kwargs)
        new = self._from_temp_dataset(ds)
        self._variable = new._variable
        self._coords = new._coords
        return self

    def compute(self, **kwargs) -> Self:
        """Trigger loading data into memory and return a new dataarray.

        Data will be computed and/or loaded from disk or a remote source.

        Unlike ``.load``, the original dataarray is left unaltered.

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
        object : DataArray
            New object with the data and all coordinates as in-memory arrays.

        See Also
        --------
        dask.compute
        DataArray.load
        DataArray.load_async
        Dataset.compute
        Variable.compute
        """
        new = self.copy(deep=False)
        return new.load(**kwargs)

    def persist(self, **kwargs) -> Self:
        """Trigger computation in constituent dask arrays

        This keeps them as dask arrays but encourages them to keep data in
        memory.  This is particularly useful when on a distributed machine.
        When on a single machine consider using ``.compute()`` instead.
        Like compute (but unlike load), the original dataset is left unaltered.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.persist``.

        Returns
        -------
        object : DataArray
            New object with all dask-backed data and coordinates as persisted dask arrays.

        See Also
        --------
        dask.persist
        """
        ds = self._to_temp_dataset().persist(**kwargs)
        return self._from_temp_dataset(ds)

    def copy(self, deep: bool = True, data: Any = None) -> Self:
        """Returns a copy of this array.

        If `deep=True`, a deep copy is made of the data array.
        Otherwise, a shallow copy is made, and the returned data array's
        values are a new view of this data array's values.

        Use `data` to create a new object with the same structure as
        original but entirely new data.

        Parameters
        ----------
        deep : bool, optional
            Whether the data array and its coordinates are loaded into memory
            and copied onto the new object. Default is True.
        data : array_like, optional
            Data to use in the new object. Must have same shape as original.
            When `data` is used, `deep` is ignored for all data variables,
            and only used for coords.

        Returns
        -------
        copy : DataArray
            New object with dimensions, attributes, coordinates, name,
            encoding, and optionally data copied from original.

        Examples
        --------
        Shallow versus deep copy

        >>> array = xr.DataArray([1, 2, 3], dims="x", coords={"x": ["a", "b", "c"]})
        >>> array.copy()
        <xarray.DataArray (x: 3)> Size: 24B
        array([1, 2, 3])
        Coordinates:
          * x        (x) <U1 12B 'a' 'b' 'c'
        >>> array_0 = array.copy(deep=False)
        >>> array_0[0] = 7
        >>> array_0
        <xarray.DataArray (x: 3)> Size: 24B
        array([7, 2, 3])
        Coordinates:
          * x        (x) <U1 12B 'a' 'b' 'c'
        >>> array
        <xarray.DataArray (x: 3)> Size: 24B
        array([7, 2, 3])
        Coordinates:
          * x        (x) <U1 12B 'a' 'b' 'c'

        Changing the data using the ``data`` argument maintains the
        structure of the original object, but with the new data. Original
        object is unaffected.

        >>> array.copy(data=[0.1, 0.2, 0.3])
        <xarray.DataArray (x: 3)> Size: 24B
        array([0.1, 0.2, 0.3])
        Coordinates:
          * x        (x) <U1 12B 'a' 'b' 'c'
        >>> array
        <xarray.DataArray (x: 3)> Size: 24B
        array([7, 2, 3])
        Coordinates:
          * x        (x) <U1 12B 'a' 'b' 'c'

        See Also
        --------
        pandas.DataFrame.copy
        """
        return self._copy(deep=deep, data=data)

    def _copy(
        self,
        deep: bool = True,
        data: Any = None,
        memo: dict[int, Any] | None = None,
    ) -> Self:
        variable = self.variable._copy(deep=deep, data=data, memo=memo)
        indexes, index_vars = self.xindexes.copy_indexes(deep=deep)

        coords = {}
        for k, v in self._coords.items():
            if k in index_vars:
                coords[k] = index_vars[k]
            else:
                coords[k] = v._copy(deep=deep, memo=memo)

        return self._replace(variable, coords, indexes=indexes)

    def __copy__(self) -> Self:
        return self._copy(deep=False)

    def __deepcopy__(self, memo: dict[int, Any] | None = None) -> Self:
        return self._copy(deep=True, memo=memo)

    # mutable objects should not be Hashable
    # https://github.com/python/mypy/issues/4266
    __hash__ = None  # type: ignore[assignment]

    @property
    def chunks(self) -> tuple[tuple[int, ...], ...] | None:
        """
        Tuple of block lengths for this dataarray's data, in order of dimensions, or None if
        the underlying data is not a dask array.

        See Also
        --------
        DataArray.chunk
        DataArray.chunksizes
        xarray.unify_chunks
        """
        return self.variable.chunks

    @property
    def chunksizes(self) -> Mapping[Any, tuple[int, ...]]:
        """
        Mapping from dimension names to block lengths for this dataarray's data.

        If this dataarray does not contain chunked arrays, the mapping will be empty.

        Cannot be modified directly, but can be modified by calling .chunk().

        Differs from DataArray.chunks because it returns a mapping of dimensions to chunk shapes
        instead of a tuple of chunk shapes.

        See Also
        --------
        DataArray.chunk
        DataArray.chunks
        xarray.unify_chunks
        """
        all_variables = [self.variable] + [c.variable for c in self.coords.values()]
        return get_chunksizes(all_variables)

    def chunk(
        self,
        chunks: T_ChunksFreq = {},  # noqa: B006  # {} even though it's technically unsafe, is being used intentionally here (#4667)
        *,
        name_prefix: str = "xarray-",
        token: str | None = None,
        lock: bool = False,
        inline_array: bool = False,
        chunked_array_type: str | ChunkManagerEntrypoint | None = None,
        from_array_kwargs=None,
        **chunks_kwargs: T_ChunkDimFreq,
    ) -> Self:
        """Coerce this array's data into a dask arrays with the given chunks.

        If this variable is a non-dask array, it will be converted to dask
        array. If it's a dask array, it will be rechunked to the given chunk
        sizes.

        If neither chunks is not provided for one or more dimensions, chunk
        sizes along that dimension will not be updated; non-dask arrays will be
        converted into dask arrays with a single block.

        Along datetime-like dimensions, a pandas frequency string is also accepted.

        Parameters
        ----------
        chunks : int, "auto", tuple of int or mapping of hashable to int or a pandas frequency string, optional
            Chunk sizes along each dimension, e.g., ``5``, ``"auto"``, ``(5, 5)`` or
            ``{"x": 5, "y": 5}`` or ``{"x": 5, "time": "YE"}``.
        name_prefix : str, optional
            Prefix for the name of the new dask array.
        token : str, optional
            Token uniquely identifying this array.
        lock : bool, default: False
            Passed on to :py:func:`dask.array.from_array`, if the array is not
            already as dask array.
        inline_array: bool, default: False
            Passed on to :py:func:`dask.array.from_array`, if the array is not
            already as dask array.
        chunked_array_type: str, optional
            Which chunked array type to coerce the underlying data array to.
            Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEntryPoint` system.
            Experimental API that should not be relied upon.
        from_array_kwargs: dict, optional
            Additional keyword arguments passed on to the `ChunkManagerEntrypoint.from_array` method used to create
            chunked arrays, via whichever chunk manager is specified through the `chunked_array_type` kwarg.
            For example, with dask as the default chunked array type, this method would pass additional kwargs
            to :py:func:`dask.array.from_array`. Experimental API that should not be relied upon.
        **chunks_kwargs : {dim: chunks, ...}, optional
            The keyword arguments form of ``chunks``.
            One of chunks or chunks_kwargs must be provided.

        Returns
        -------
        chunked : xarray.DataArray

        See Also
        --------
        DataArray.chunks
        DataArray.chunksizes
        xarray.unify_chunks
        dask.array.from_array
        """
        chunk_mapping: T_ChunksFreq
        if chunks is None:
            warnings.warn(
                "None value for 'chunks' is deprecated. "
                "It will raise an error in the future. Use instead '{}'",
                category=FutureWarning,
                stacklevel=2,
            )
            chunk_mapping = {}

        if isinstance(chunks, float | str | int):
            # ignoring type; unclear why it won't accept a Literal into the value.
            chunk_mapping = dict.fromkeys(self.dims, chunks)
        elif isinstance(chunks, tuple | list):
            utils.emit_user_level_warning(
                "Supplying chunks as dimension-order tuples is deprecated. "
                "It will raise an error in the future. Instead use a dict with dimension names as keys.",
                category=DeprecationWarning,
            )
            if len(chunks) != len(self.dims):
                raise ValueError(
                    f"chunks must have the same number of elements as dimensions. "
                    f"Expected {len(self.dims)} elements, got {len(chunks)}."
                )
            chunk_mapping = dict(zip(self.dims, chunks, strict=True))
        else:
            chunk_mapping = either_dict_or_kwargs(chunks, chunks_kwargs, "chunk")

        ds = self._to_temp_dataset().chunk(
            chunk_mapping,
            name_prefix=name_prefix,
            token=token,
            lock=lock,
            inline_array=inline_array,
            chunked_array_type=chunked_array_type,
            from_array_kwargs=from_array_kwargs,
        )
        return self._from_temp_dataset(ds)

    def isel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        drop: bool = False,
        missing_dims: ErrorOptionsWithWarn = "raise",
        **indexers_kwargs: Any,
    ) -> Self:
        """Return a new DataArray whose data is given by selecting indexes
        along the specified dimension(s).

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
            DataArray:
            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions
        **indexers_kwargs : {dim: indexer, ...}, optional
            The keyword arguments form of ``indexers``.

        Returns
        -------
        indexed : xarray.DataArray

        See Also
        --------
        :func:`Dataset.isel <Dataset.isel>`
        :func:`DataArray.sel <DataArray.sel>`

        :doc:`xarray-tutorial:intermediate/indexing/indexing`
            Tutorial material on indexing with Xarray objects

        :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
            Tutorial material on basics of indexing

        Examples
        --------
        >>> da = xr.DataArray(np.arange(25).reshape(5, 5), dims=("x", "y"))
        >>> da
        <xarray.DataArray (x: 5, y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Dimensions without coordinates: x, y

        >>> tgt_x = xr.DataArray(np.arange(0, 5), dims="points")
        >>> tgt_y = xr.DataArray(np.arange(0, 5), dims="points")
        >>> da = da.isel(x=tgt_x, y=tgt_y)
        >>> da
        <xarray.DataArray (points: 5)> Size: 40B
        array([ 0,  6, 12, 18, 24])
        Dimensions without coordinates: points
        """

        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "isel")

        if any(is_fancy_indexer(idx) for idx in indexers.values()):
            ds = self._to_temp_dataset()._isel_fancy(
                indexers, drop=drop, missing_dims=missing_dims
            )
            return self._from_temp_dataset(ds)

        # Much faster algorithm for when all indexers are ints, slices, one-dimensional
        # lists, or zero or one-dimensional np.ndarray's

        variable = self._variable.isel(indexers, missing_dims=missing_dims)
        indexes, index_variables = isel_indexes(self.xindexes, indexers)

        coords = {}
        for coord_name, coord_value in self._coords.items():
            if coord_name in index_variables:
                coord_value = index_variables[coord_name]
            else:
                coord_indexers = {
                    k: v for k, v in indexers.items() if k in coord_value.dims
                }
                if coord_indexers:
                    coord_value = coord_value.isel(coord_indexers)
                    if drop and coord_value.ndim == 0:
                        continue
            coords[coord_name] = coord_value

        return self._replace(variable=variable, coords=coords, indexes=indexes)

    def sel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        method: str | None = None,
        tolerance=None,
        drop: bool = False,
        **indexers_kwargs: Any,
    ) -> Self:
        """Return a new DataArray whose data is given by selecting index
        labels along the specified dimension(s).

        In contrast to `DataArray.isel`, indexers for this method should use
        labels instead of integers.

        Under the hood, this method is powered by using pandas's powerful Index
        objects. This makes label based indexing essentially just as fast as
        using integer indexing.

        It also means this method uses pandas's (well documented) logic for
        indexing. This means you can use string shortcuts for datetime indexes
        (e.g., '2000-01' to select all values in January 2000). It also means
        that slices are treated as inclusive of both the start and stop values,
        unlike normal Python indexing.

        .. warning::

          Do not try to assign values when using any of the indexing methods
          ``isel`` or ``sel``::

            da = xr.DataArray([0, 1, 2, 3], dims=["x"])
            # DO NOT do this
            da.isel(x=[0, 1, 2])[1] = -1

          Assigning values with the chained indexing using ``.sel`` or
          ``.isel`` fails silently.

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

            - None (default): only exact matches
            - pad / ffill: propagate last valid index value forward
            - backfill / bfill: propagate next valid index value backward
            - nearest: use nearest valid index value

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
        obj : DataArray
            A new DataArray with the same contents as this DataArray, except the
            data and each dimension is indexed by the appropriate indexers.
            If indexer DataArrays have coordinates that do not conflict with
            this object, then these coordinates will be attached.
            In general, each array's data will be a view of the array's data
            in this DataArray, unless vectorized indexing was triggered by using
            an array indexer, in which case the data will be a copy.

        See Also
        --------
        :func:`Dataset.sel <Dataset.sel>`
        :func:`DataArray.isel <DataArray.isel>`

        :doc:`xarray-tutorial:intermediate/indexing/indexing`
            Tutorial material on indexing with Xarray objects

        :doc:`xarray-tutorial:fundamentals/02.1_indexing_Basic`
            Tutorial material on basics of indexing

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.arange(25).reshape(5, 5),
        ...     coords={"x": np.arange(5), "y": np.arange(5)},
        ...     dims=("x", "y"),
        ... )
        >>> da
        <xarray.DataArray (x: 5, y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Coordinates:
          * x        (x) int64 40B 0 1 2 3 4
          * y        (y) int64 40B 0 1 2 3 4

        >>> tgt_x = xr.DataArray(np.linspace(0, 4, num=5), dims="points")
        >>> tgt_y = xr.DataArray(np.linspace(0, 4, num=5), dims="points")
        >>> da = da.sel(x=tgt_x, y=tgt_y, method="nearest")
        >>> da
        <xarray.DataArray (points: 5)> Size: 40B
        array([ 0,  6, 12, 18, 24])
        Coordinates:
            x        (points) int64 40B 0 1 2 3 4
            y        (points) int64 40B 0 1 2 3 4
        Dimensions without coordinates: points
        """
        ds = self._to_temp_dataset().sel(
            indexers=indexers,
            drop=drop,
            method=method,
            tolerance=tolerance,
            **indexers_kwargs,
        )
        return self._from_temp_dataset(ds)

    def _shuffle(
        self, dim: Hashable, *, indices: GroupIndices, chunks: T_Chunks
    ) -> Self:
        ds = self._to_temp_dataset()._shuffle(dim=dim, indices=indices, chunks=chunks)
        return self._from_temp_dataset(ds)

    def head(
        self,
        indexers: Mapping[Any, int] | int | None = None,
        **indexers_kwargs: Any,
    ) -> Self:
        """Return a new DataArray whose data is given by the the first `n`
        values along the specified dimension(s). Default `n` = 5

        See Also
        --------
        Dataset.head
        DataArray.tail
        DataArray.thin

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.arange(25).reshape(5, 5),
        ...     dims=("x", "y"),
        ... )
        >>> da
        <xarray.DataArray (x: 5, y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Dimensions without coordinates: x, y

        >>> da.head(x=1)
        <xarray.DataArray (x: 1, y: 5)> Size: 40B
        array([[0, 1, 2, 3, 4]])
        Dimensions without coordinates: x, y

        >>> da.head({"x": 2, "y": 2})
        <xarray.DataArray (x: 2, y: 2)> Size: 32B
        array([[0, 1],
               [5, 6]])
        Dimensions without coordinates: x, y
        """
        ds = self._to_temp_dataset().head(indexers, **indexers_kwargs)
        return self._from_temp_dataset(ds)

    def tail(
        self,
        indexers: Mapping[Any, int] | int | None = None,
        **indexers_kwargs: Any,
    ) -> Self:
        """Return a new DataArray whose data is given by the the last `n`
        values along the specified dimension(s). Default `n` = 5

        See Also
        --------
        Dataset.tail
        DataArray.head
        DataArray.thin

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.arange(25).reshape(5, 5),
        ...     dims=("x", "y"),
        ... )
        >>> da
        <xarray.DataArray (x: 5, y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Dimensions without coordinates: x, y

        >>> da.tail(y=1)
        <xarray.DataArray (x: 5, y: 1)> Size: 40B
        array([[ 4],
               [ 9],
               [14],
               [19],
               [24]])
        Dimensions without coordinates: x, y

        >>> da.tail({"x": 2, "y": 2})
        <xarray.DataArray (x: 2, y: 2)> Size: 32B
        array([[18, 19],
               [23, 24]])
        Dimensions without coordinates: x, y
        """
        ds = self._to_temp_dataset().tail(indexers, **indexers_kwargs)
        return self._from_temp_dataset(ds)

    def thin(
        self,
        indexers: Mapping[Any, int] | int | None = None,
        **indexers_kwargs: Any,
    ) -> Self:
        """Return a new DataArray whose data is given by each `n` value
        along the specified dimension(s).

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
        >>> x
        <xarray.DataArray (x: 2, y: 13)> Size: 208B
        array([[ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12],
               [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]])
        Coordinates:
          * x        (x) int64 16B 0 1
          * y        (y) int64 104B 0 1 2 3 4 5 6 7 8 9 10 11 12

        >>>
        >>> x.thin(3)
        <xarray.DataArray (x: 1, y: 5)> Size: 40B
        array([[ 0,  3,  6,  9, 12]])
        Coordinates:
          * x        (x) int64 8B 0
          * y        (y) int64 40B 0 3 6 9 12
        >>> x.thin({"x": 2, "y": 5})
        <xarray.DataArray (x: 1, y: 3)> Size: 24B
        array([[ 0,  5, 10]])
        Coordinates:
          * x        (x) int64 8B 0
          * y        (y) int64 24B 0 5 10

        See Also
        --------
        Dataset.thin
        DataArray.head
        DataArray.tail
        """
        ds = self._to_temp_dataset().thin(indexers, **indexers_kwargs)
        return self._from_temp_dataset(ds)

    def broadcast_like(
        self,
        other: T_DataArrayOrSet,
        *,
        exclude: Iterable[Hashable] | None = None,
    ) -> Self:
        """Broadcast this DataArray against another Dataset or DataArray.

        This is equivalent to xr.broadcast(other, self)[1]

        xarray objects are broadcast against each other in arithmetic
        operations, so this method is not be necessary for most uses.

        If no change is needed, the input data is returned to the output
        without being copied.

        If new coords are added by the broadcast, their values are
        NaN filled.

        Parameters
        ----------
        other : Dataset or DataArray
            Object against which to broadcast this array.
        exclude : iterable of Hashable, optional
            Dimensions that must not be broadcasted

        Returns
        -------
        new_da : DataArray
            The caller broadcasted against ``other``.

        Examples
        --------
        >>> arr1 = xr.DataArray(
        ...     np.random.randn(2, 3),
        ...     dims=("x", "y"),
        ...     coords={"x": ["a", "b"], "y": ["a", "b", "c"]},
        ... )
        >>> arr2 = xr.DataArray(
        ...     np.random.randn(3, 2),
        ...     dims=("x", "y"),
        ...     coords={"x": ["a", "b", "c"], "y": ["a", "b"]},
        ... )
        >>> arr1
        <xarray.DataArray (x: 2, y: 3)> Size: 48B
        array([[ 1.76405235,  0.40015721,  0.97873798],
               [ 2.2408932 ,  1.86755799, -0.97727788]])
        Coordinates:
          * x        (x) <U1 8B 'a' 'b'
          * y        (y) <U1 12B 'a' 'b' 'c'
        >>> arr2
        <xarray.DataArray (x: 3, y: 2)> Size: 48B
        array([[ 0.95008842, -0.15135721],
               [-0.10321885,  0.4105985 ],
               [ 0.14404357,  1.45427351]])
        Coordinates:
          * x        (x) <U1 12B 'a' 'b' 'c'
          * y        (y) <U1 8B 'a' 'b'
        >>> arr1.broadcast_like(arr2)
        <xarray.DataArray (x: 3, y: 3)> Size: 72B
        array([[ 1.76405235,  0.40015721,  0.97873798],
               [ 2.2408932 ,  1.86755799, -0.97727788],
               [        nan,         nan,         nan]])
        Coordinates:
          * x        (x) <U1 12B 'a' 'b' 'c'
          * y        (y) <U1 12B 'a' 'b' 'c'
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
        """Callback called from ``Aligner`` to create a new reindexed DataArray."""

        if isinstance(fill_value, dict):
            fill_value = fill_value.copy()
            sentinel = object()
            value = fill_value.pop(self.name, sentinel)
            if value is not sentinel:
                fill_value[_THIS_ARRAY] = value

        ds = self._to_temp_dataset()
        reindexed = ds._reindex_callback(
            aligner,
            dim_pos_indexers,
            variables,
            indexes,
            fill_value,
            exclude_dims,
            exclude_vars,
        )

        da = self._from_temp_dataset(reindexed)
        da.encoding = self.encoding

        return da

    def reindex_like(
        self,
        other: T_DataArrayOrSet,
        *,
        method: ReindexMethodOptions = None,
        tolerance: float | Iterable[float] | str | None = None,
        copy: bool = True,
        fill_value=dtypes.NA,
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
        method : {None, "nearest", "pad", "ffill", "backfill", "bfill"}, optional
            Method to use for filling index values from other not found on this
            data array:

            - None (default): don't fill gaps
            - pad / ffill: propagate last valid index value forward
            - backfill / bfill: propagate next valid index value backward
            - nearest: use nearest valid index value

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
            Value to use for newly missing values. If a dict-like, maps
            variable names (including coordinates) to fill values. Use this
            data array's name to refer to the data array's values.

        Returns
        -------
        reindexed : DataArray
            Another dataset array, with this array's data but coordinates from
            the other object.

        Examples
        --------
        >>> data = np.arange(12).reshape(4, 3)
        >>> da1 = xr.DataArray(
        ...     data=data,
        ...     dims=["x", "y"],
        ...     coords={"x": [10, 20, 30, 40], "y": [70, 80, 90]},
        ... )
        >>> da1
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * x        (x) int64 32B 10 20 30 40
          * y        (y) int64 24B 70 80 90
        >>> da2 = xr.DataArray(
        ...     data=data,
        ...     dims=["x", "y"],
        ...     coords={"x": [40, 30, 20, 10], "y": [90, 80, 70]},
        ... )
        >>> da2
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * x        (x) int64 32B 40 30 20 10
          * y        (y) int64 24B 90 80 70

        Reindexing with both DataArrays having the same coordinates set, but in different order:

        >>> da1.reindex_like(da2)
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[11, 10,  9],
               [ 8,  7,  6],
               [ 5,  4,  3],
               [ 2,  1,  0]])
        Coordinates:
          * x        (x) int64 32B 40 30 20 10
          * y        (y) int64 24B 90 80 70

        Reindexing with the other array having additional coordinates:

        >>> da3 = xr.DataArray(
        ...     data=data,
        ...     dims=["x", "y"],
        ...     coords={"x": [20, 10, 29, 39], "y": [70, 80, 90]},
        ... )
        >>> da1.reindex_like(da3)
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 3.,  4.,  5.],
               [ 0.,  1.,  2.],
               [nan, nan, nan],
               [nan, nan, nan]])
        Coordinates:
          * x        (x) int64 32B 20 10 29 39
          * y        (y) int64 24B 70 80 90

        Filling missing values with the previous valid index with respect to the coordinates' value:

        >>> da1.reindex_like(da3, method="ffill")
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[3, 4, 5],
               [0, 1, 2],
               [3, 4, 5],
               [6, 7, 8]])
        Coordinates:
          * x        (x) int64 32B 20 10 29 39
          * y        (y) int64 24B 70 80 90

        Filling missing values while tolerating specified error for inexact matches:

        >>> da1.reindex_like(da3, method="ffill", tolerance=5)
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 3.,  4.,  5.],
               [ 0.,  1.,  2.],
               [nan, nan, nan],
               [nan, nan, nan]])
        Coordinates:
          * x        (x) int64 32B 20 10 29 39
          * y        (y) int64 24B 70 80 90

        Filling missing values with manually specified values:

        >>> da1.reindex_like(da3, fill_value=19)
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 3,  4,  5],
               [ 0,  1,  2],
               [19, 19, 19],
               [19, 19, 19]])
        Coordinates:
          * x        (x) int64 32B 20 10 29 39
          * y        (y) int64 24B 70 80 90

        Note that unlike ``broadcast_like``, ``reindex_like`` doesn't create new dimensions:

        >>> da1.sel(x=20)
        <xarray.DataArray (y: 3)> Size: 24B
        array([3, 4, 5])
        Coordinates:
            x        int64 8B 20
          * y        (y) int64 24B 70 80 90

        ...so ``b`` in not added here:

        >>> da1.sel(x=20).reindex_like(da1)
        <xarray.DataArray (y: 3)> Size: 24B
        array([3, 4, 5])
        Coordinates:
            x        int64 8B 20
          * y        (y) int64 24B 70 80 90

        See Also
        --------
        DataArray.reindex
        DataArray.broadcast_like
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
        *,
        method: ReindexMethodOptions = None,
        tolerance: float | Iterable[float] | str | None = None,
        copy: bool = True,
        fill_value=dtypes.NA,
        **indexers_kwargs: Any,
    ) -> Self:
        """Conform this object onto the indexes of another object, filling in
        missing values with ``fill_value``. The default fill value is NaN.

        Parameters
        ----------
        indexers : dict, optional
            Dictionary with keys given by dimension names and values given by
            arrays of coordinates tick labels. Any mismatched coordinate
            values will be filled in with NaN, and any mismatched dimension
            names will simply be ignored.
            One of indexers or indexers_kwargs must be provided.
        copy : bool, optional
            If ``copy=True``, data in the return value is always copied. If
            ``copy=False`` and reindexing is unnecessary, or can be performed
            with only slice operations, then the output may share memory with
            the input. In either case, a new xarray object is always returned.
        method : {None, 'nearest', 'pad'/'ffill', 'backfill'/'bfill'}, optional
            Method to use for filling index values in ``indexers`` not found on
            this data array:

            - None (default): don't fill gaps
            - pad / ffill: propagate last valid index value forward
            - backfill / bfill: propagate next valid index value backward
            - nearest: use nearest valid index value

        tolerance : float | Iterable[float] | str | None, default: None
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations must
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.
            Tolerance may be a scalar value, which applies the same tolerance
            to all values, or list-like, which applies variable tolerance per
            element. List-like must be the same size as the index and its dtype
            must exactly match the index’s type.
        fill_value : scalar or dict-like, optional
            Value to use for newly missing values. If a dict-like, maps
            variable names (including coordinates) to fill values. Use this
            data array's name to refer to the data array's values.
        **indexers_kwargs : {dim: indexer, ...}, optional
            The keyword arguments form of ``indexers``.
            One of indexers or indexers_kwargs must be provided.

        Returns
        -------
        reindexed : DataArray
            Another dataset array, with this array's data but replaced
            coordinates.

        Examples
        --------
        Reverse latitude:

        >>> da = xr.DataArray(
        ...     np.arange(4),
        ...     coords=[np.array([90, 89, 88, 87])],
        ...     dims="lat",
        ... )
        >>> da
        <xarray.DataArray (lat: 4)> Size: 32B
        array([0, 1, 2, 3])
        Coordinates:
          * lat      (lat) int64 32B 90 89 88 87
        >>> da.reindex(lat=da.lat[::-1])
        <xarray.DataArray (lat: 4)> Size: 32B
        array([3, 2, 1, 0])
        Coordinates:
          * lat      (lat) int64 32B 87 88 89 90

        See Also
        --------
        DataArray.reindex_like
        align
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

    def interp(
        self,
        coords: Mapping[Any, Any] | None = None,
        method: InterpOptions = "linear",
        assume_sorted: bool = False,
        kwargs: Mapping[str, Any] | None = None,
        **coords_kwargs: Any,
    ) -> Self:
        """
        Interpolate a DataArray onto new coordinates.

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
            If False, values of x can be in any order and they are sorted
            first. If True, x has to be an array of monotonically increasing
            values.
        kwargs : dict-like or None, default: None
            Additional keyword arguments passed to scipy's interpolator. Valid
            options and their behavior depend whether ``interp1d`` or
            ``interpn`` is used.
        **coords_kwargs : {dim: coordinate, ...}, optional
            The keyword arguments form of ``coords``.
            One of coords or coords_kwargs must be provided.

        Returns
        -------
        interpolated : DataArray
            New dataarray on the new coordinates.

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
            Tutorial material on manipulating data resolution using :py:func:`~xarray.DataArray.interp`

        Examples
        --------
        >>> da = xr.DataArray(
        ...     data=[[1, 4, 2, 9], [2, 7, 6, np.nan], [6, np.nan, 5, 8]],
        ...     dims=("x", "y"),
        ...     coords={"x": [0, 1, 2], "y": [10, 12, 14, 16]},
        ... )
        >>> da
        <xarray.DataArray (x: 3, y: 4)> Size: 96B
        array([[ 1.,  4.,  2.,  9.],
               [ 2.,  7.,  6., nan],
               [ 6., nan,  5.,  8.]])
        Coordinates:
          * x        (x) int64 24B 0 1 2
          * y        (y) int64 32B 10 12 14 16

        1D linear interpolation (the default):

        >>> da.interp(x=[0, 0.75, 1.25, 1.75])
        <xarray.DataArray (x: 4, y: 4)> Size: 128B
        array([[1.  , 4.  , 2.  ,  nan],
               [1.75, 6.25, 5.  ,  nan],
               [3.  ,  nan, 5.75,  nan],
               [5.  ,  nan, 5.25,  nan]])
        Coordinates:
          * y        (y) int64 32B 10 12 14 16
          * x        (x) float64 32B 0.0 0.75 1.25 1.75

        1D nearest interpolation:

        >>> da.interp(x=[0, 0.75, 1.25, 1.75], method="nearest")
        <xarray.DataArray (x: 4, y: 4)> Size: 128B
        array([[ 1.,  4.,  2.,  9.],
               [ 2.,  7.,  6., nan],
               [ 2.,  7.,  6., nan],
               [ 6., nan,  5.,  8.]])
        Coordinates:
          * y        (y) int64 32B 10 12 14 16
          * x        (x) float64 32B 0.0 0.75 1.25 1.75

        1D linear extrapolation:

        >>> da.interp(
        ...     x=[1, 1.5, 2.5, 3.5],
        ...     method="linear",
        ...     kwargs={"fill_value": "extrapolate"},
        ... )
        <xarray.DataArray (x: 4, y: 4)> Size: 128B
        array([[ 2. ,  7. ,  6. ,  nan],
               [ 4. ,  nan,  5.5,  nan],
               [ 8. ,  nan,  4.5,  nan],
               [12. ,  nan,  3.5,  nan]])
        Coordinates:
          * y        (y) int64 32B 10 12 14 16
          * x        (x) float64 32B 1.0 1.5 2.5 3.5

        2D linear interpolation:

        >>> da.interp(x=[0, 0.75, 1.25, 1.75], y=[11, 13, 15], method="linear")
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[2.5  , 3.   ,   nan],
               [4.   , 5.625,   nan],
               [  nan,   nan,   nan],
               [  nan,   nan,   nan]])
        Coordinates:
          * x        (x) float64 32B 0.0 0.75 1.25 1.75
          * y        (y) int64 24B 11 13 15
        """
        if self.dtype.kind not in "uifc":
            raise TypeError(
                f"interp only works for a numeric type array. Given {self.dtype}."
            )
        ds = self._to_temp_dataset().interp(
            coords,
            method=method,
            kwargs=kwargs,
            assume_sorted=assume_sorted,
            **coords_kwargs,
        )
        return self._from_temp_dataset(ds)

    def interp_like(
        self,
        other: T_Xarray,
        method: InterpOptions = "linear",
        assume_sorted: bool = False,
        kwargs: Mapping[str, Any] | None = None,
    ) -> Self:
        """Interpolate this object onto the coordinates of another object,
        filling out of range values with NaN.

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
            Additional keyword arguments passed to the interpolant.

        Returns
        -------
        interpolated : DataArray
            Another dataarray by interpolating this dataarray's data along the
            coordinates of the other object.

        Notes
        -----
        - scipy is required.
        - If the dataarray has object-type coordinates, reindex is used for these
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
        :func:`DataArray.interp`
        :func:`DataArray.reindex_like`
        :mod:`scipy.interpolate`

        Examples
        --------
        >>> data = np.arange(12).reshape(4, 3)
        >>> da1 = xr.DataArray(
        ...     data=data,
        ...     dims=["x", "y"],
        ...     coords={"x": [10, 20, 30, 40], "y": [70, 80, 90]},
        ... )
        >>> da1
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * x        (x) int64 32B 10 20 30 40
          * y        (y) int64 24B 70 80 90
        >>> da2 = xr.DataArray(
        ...     data=data,
        ...     dims=["x", "y"],
        ...     coords={"x": [10, 20, 29, 39], "y": [70, 80, 90]},
        ... )
        >>> da2
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * x        (x) int64 32B 10 20 29 39
          * y        (y) int64 24B 70 80 90

        Interpolate the values in the coordinates of the other DataArray with respect to the source's values:

        >>> da2.interp_like(da1)
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[0. , 1. , 2. ],
               [3. , 4. , 5. ],
               [6.3, 7.3, 8.3],
               [nan, nan, nan]])
        Coordinates:
          * x        (x) int64 32B 10 20 30 40
          * y        (y) int64 24B 70 80 90

        Could also extrapolate missing values:

        >>> da2.interp_like(da1, kwargs={"fill_value": "extrapolate"})
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0. ,  1. ,  2. ],
               [ 3. ,  4. ,  5. ],
               [ 6.3,  7.3,  8.3],
               [ 9.3, 10.3, 11.3]])
        Coordinates:
          * x        (x) int64 32B 10 20 30 40
          * y        (y) int64 24B 70 80 90
        """

        if self.dtype.kind not in "uifc":
            raise TypeError(
                f"interp only works for a numeric type array. Given {self.dtype}."
            )
        ds = self._to_temp_dataset().interp_like(
            other, method=method, kwargs=kwargs, assume_sorted=assume_sorted
        )
        return self._from_temp_dataset(ds)

    def rename(
        self,
        new_name_or_name_dict: Hashable | Mapping[Any, Hashable] | None = None,
        **names: Hashable,
    ) -> Self:
        """Returns a new DataArray with renamed coordinates, dimensions or a new name.

        Parameters
        ----------
        new_name_or_name_dict : str or dict-like, optional
            If the argument is dict-like, it used as a mapping from old
            names to new names for coordinates or dimensions. Otherwise,
            use the argument as the new name for this array.
        **names : Hashable, optional
            The keyword arguments form of a mapping from old names to
            new names for coordinates or dimensions.
            One of new_name_or_name_dict or names must be provided.

        Returns
        -------
        renamed : DataArray
            Renamed array or array with renamed coordinates.

        See Also
        --------
        Dataset.rename
        DataArray.swap_dims
        """
        if new_name_or_name_dict is None and not names:
            # change name to None?
            return self._replace(name=None)
        if utils.is_dict_like(new_name_or_name_dict) or new_name_or_name_dict is None:
            # change dims/coords
            name_dict = either_dict_or_kwargs(new_name_or_name_dict, names, "rename")
            dataset = self._to_temp_dataset()._rename(name_dict)
            return self._from_temp_dataset(dataset)
        if utils.hashable(new_name_or_name_dict) and names:
            # change name + dims/coords
            dataset = self._to_temp_dataset()._rename(names)
            dataarray = self._from_temp_dataset(dataset)
            return dataarray._replace(name=new_name_or_name_dict)
        # only change name
        return self._replace(name=new_name_or_name_dict)

    def swap_dims(
        self,
        dims_dict: Mapping[Any, Hashable] | None = None,
        **dims_kwargs,
    ) -> Self:
        """Returns a new DataArray with swapped dimensions.

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
        swapped : DataArray
            DataArray with swapped dimensions.

        Examples
        --------
        >>> arr = xr.DataArray(
        ...     data=[0, 1],
        ...     dims="x",
        ...     coords={"x": ["a", "b"], "y": ("x", [0, 1])},
        ... )
        >>> arr
        <xarray.DataArray (x: 2)> Size: 16B
        array([0, 1])
        Coordinates:
          * x        (x) <U1 8B 'a' 'b'
            y        (x) int64 16B 0 1

        >>> arr.swap_dims({"x": "y"})
        <xarray.DataArray (y: 2)> Size: 16B
        array([0, 1])
        Coordinates:
            x        (y) <U1 8B 'a' 'b'
          * y        (y) int64 16B 0 1

        >>> arr.swap_dims({"x": "z"})
        <xarray.DataArray (z: 2)> Size: 16B
        array([0, 1])
        Coordinates:
            x        (z) <U1 8B 'a' 'b'
            y        (z) int64 16B 0 1
        Dimensions without coordinates: z

        See Also
        --------
        DataArray.rename
        Dataset.swap_dims
        """
        dims_dict = either_dict_or_kwargs(dims_dict, dims_kwargs, "swap_dims")
        ds = self._to_temp_dataset().swap_dims(dims_dict)
        return self._from_temp_dataset(ds)

    def expand_dims(
        self,
        dim: Hashable | Sequence[Hashable] | Mapping[Any, Any] | None = None,
        axis: int | Sequence[int] | None = None,
        create_index_for_new_dim: bool = True,
        **dim_kwargs: Any,
    ) -> Self:
        """Return a new object with an additional axis (or axes) inserted at
        the corresponding position in the array shape. The new object is a
        view into the underlying array, not a copy.

        If dim is already a scalar coordinate, it will be promoted to a 1D
        coordinate consisting of a single value.

        The automatic creation of indexes to back new 1D coordinate variables
        controlled by the create_index_for_new_dim kwarg.

        Parameters
        ----------
        dim : Hashable, sequence of Hashable, dict, or None, optional
            Dimensions to include on the new variable.
            If provided as str or sequence of str, then dimensions are inserted
            with length 1. If provided as a dict, then the keys are the new
            dimensions and the values are either integers (giving the length of
            the new dimensions) or sequence/ndarray (giving the coordinates of
            the new dimensions).
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
        expanded : DataArray
            This object, but with additional dimension(s).

        See Also
        --------
        Dataset.expand_dims

        Examples
        --------
        >>> da = xr.DataArray(np.arange(5), dims=("x"))
        >>> da
        <xarray.DataArray (x: 5)> Size: 40B
        array([0, 1, 2, 3, 4])
        Dimensions without coordinates: x

        Add new dimension of length 2:

        >>> da.expand_dims(dim={"y": 2})
        <xarray.DataArray (y: 2, x: 5)> Size: 80B
        array([[0, 1, 2, 3, 4],
               [0, 1, 2, 3, 4]])
        Dimensions without coordinates: y, x

        >>> da.expand_dims(dim={"y": 2}, axis=1)
        <xarray.DataArray (x: 5, y: 2)> Size: 80B
        array([[0, 0],
               [1, 1],
               [2, 2],
               [3, 3],
               [4, 4]])
        Dimensions without coordinates: x, y

        Add a new dimension with coordinates from array:

        >>> da.expand_dims(dim={"y": np.arange(5)}, axis=0)
        <xarray.DataArray (y: 5, x: 5)> Size: 200B
        array([[0, 1, 2, 3, 4],
               [0, 1, 2, 3, 4],
               [0, 1, 2, 3, 4],
               [0, 1, 2, 3, 4],
               [0, 1, 2, 3, 4]])
        Coordinates:
          * y        (y) int64 40B 0 1 2 3 4
        Dimensions without coordinates: x
        """
        if isinstance(dim, int):
            raise TypeError("dim should be Hashable or sequence/mapping of Hashables")
        elif isinstance(dim, Sequence) and not isinstance(dim, str):
            if len(dim) != len(set(dim)):
                raise ValueError("dims should not contain duplicate values.")
            dim = dict.fromkeys(dim, 1)
        elif dim is not None and not isinstance(dim, Mapping):
            dim = {dim: 1}

        dim = either_dict_or_kwargs(dim, dim_kwargs, "expand_dims")
        ds = self._to_temp_dataset().expand_dims(
            dim, axis, create_index_for_new_dim=create_index_for_new_dim
        )
        return self._from_temp_dataset(ds)

    def set_index(
        self,
        indexes: Mapping[Any, Hashable | Sequence[Hashable]] | None = None,
        append: bool = False,
        **indexes_kwargs: Hashable | Sequence[Hashable],
    ) -> Self:
        """Set DataArray (multi-)indexes using one or more existing
        coordinates.

        This legacy method is limited to pandas (multi-)indexes and
        1-dimensional "dimension" coordinates. See
        :py:meth:`~DataArray.set_xindex` for setting a pandas or a custom
        Xarray-compatible index from one or more arbitrary coordinates.

        Parameters
        ----------
        indexes : {dim: index, ...}
            Mapping from names matching dimensions and values given
            by (lists of) the names of existing coordinates or variables to set
            as new (multi-)index.
        append : bool, default: False
            If True, append the supplied index(es) to the existing index(es).
            Otherwise replace the existing index(es).
        **indexes_kwargs : optional
            The keyword arguments form of ``indexes``.
            One of indexes or indexes_kwargs must be provided.

        Returns
        -------
        obj : DataArray
            Another DataArray, with this data but replaced coordinates.

        Examples
        --------
        >>> arr = xr.DataArray(
        ...     data=np.ones((2, 3)),
        ...     dims=["x", "y"],
        ...     coords={"x": range(2), "y": range(3), "a": ("x", [3, 4])},
        ... )
        >>> arr
        <xarray.DataArray (x: 2, y: 3)> Size: 48B
        array([[1., 1., 1.],
               [1., 1., 1.]])
        Coordinates:
          * x        (x) int64 16B 0 1
          * y        (y) int64 24B 0 1 2
            a        (x) int64 16B 3 4
        >>> arr.set_index(x="a")
        <xarray.DataArray (x: 2, y: 3)> Size: 48B
        array([[1., 1., 1.],
               [1., 1., 1.]])
        Coordinates:
          * x        (x) int64 16B 3 4
          * y        (y) int64 24B 0 1 2

        See Also
        --------
        DataArray.reset_index
        DataArray.set_xindex
        """
        ds = self._to_temp_dataset().set_index(indexes, append=append, **indexes_kwargs)
        return self._from_temp_dataset(ds)

    def reset_index(
        self,
        dims_or_levels: Hashable | Sequence[Hashable],
        drop: bool = False,
    ) -> Self:
        """Reset the specified index(es) or multi-index level(s).

        This legacy method is specific to pandas (multi-)indexes and
        1-dimensional "dimension" coordinates. See the more generic
        :py:meth:`~DataArray.drop_indexes` and :py:meth:`~DataArray.set_xindex`
        method to respectively drop and set pandas or custom indexes for
        arbitrary coordinates.

        Parameters
        ----------
        dims_or_levels : Hashable or sequence of Hashable
            Name(s) of the dimension(s) and/or multi-index level(s) that will
            be reset.
        drop : bool, default: False
            If True, remove the specified indexes and/or multi-index levels
            instead of extracting them as new coordinates (default: False).

        Returns
        -------
        obj : DataArray
            Another dataarray, with this dataarray's data but replaced
            coordinates.

        See Also
        --------
        DataArray.set_index
        DataArray.set_xindex
        DataArray.drop_indexes
        """
        ds = self._to_temp_dataset().reset_index(dims_or_levels, drop=drop)
        return self._from_temp_dataset(ds)

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
        index_cls : subclass of :class:`~xarray.indexes.Index`
            The type of index to create. By default, try setting
            a pandas (multi-)index from the supplied coordinates.
        **options
            Options passed to the index constructor.

        Returns
        -------
        obj : DataArray
            Another dataarray, with this dataarray's data and with a new index.

        """
        ds = self._to_temp_dataset().set_xindex(coord_names, index_cls, **options)
        return self._from_temp_dataset(ds)

    def reorder_levels(
        self,
        dim_order: Mapping[Any, Sequence[int | Hashable]] | None = None,
        **dim_order_kwargs: Sequence[int | Hashable],
    ) -> Self:
        """Rearrange index levels using input order.

        Parameters
        ----------
        dim_order dict-like of Hashable to int or Hashable: optional
            Mapping from names matching dimensions and values given
            by lists representing new level orders. Every given dimension
            must have a multi-index.
        **dim_order_kwargs : optional
            The keyword arguments form of ``dim_order``.
            One of dim_order or dim_order_kwargs must be provided.

        Returns
        -------
        obj : DataArray
            Another dataarray, with this dataarray's data but replaced
            coordinates.
        """
        ds = self._to_temp_dataset().reorder_levels(dim_order, **dim_order_kwargs)
        return self._from_temp_dataset(ds)

    @partial(deprecate_dims, old_name="dimensions")
    def stack(
        self,
        dim: Mapping[Any, Sequence[Hashable]] | None = None,
        create_index: bool | None = True,
        index_cls: type[Index] = PandasMultiIndex,
        **dim_kwargs: Sequence[Hashable | EllipsisType],
    ) -> Self:
        """
        Stack any number of existing dimensions into a single new dimension.

        New dimensions will be added at the end, and the corresponding
        coordinate variables will be combined into a MultiIndex.

        Parameters
        ----------
        dim : mapping of Hashable to sequence of Hashable
            Mapping of the form `new_name=(dim1, dim2, ...)`.
            Names of new dimensions, and the existing dimensions that they
            replace. An ellipsis (`...`) will be replaced by all unlisted dimensions.
            Passing a list containing an ellipsis (`stacked_dim=[...]`) will stack over
            all dimensions.
        create_index : bool or None, default: True
            If True, create a multi-index for each of the stacked dimensions.
            If False, don't create any index.
            If None, create a multi-index only if exactly one single (1-d) coordinate
            index is found for every dimension to stack.
        index_cls: class, optional
            Can be used to pass a custom multi-index type. Must be an Xarray index that
            implements `.stack()`. By default, a pandas multi-index wrapper is used.
        **dim_kwargs
            The keyword arguments form of ``dim``.
            One of dim or dim_kwargs must be provided.

        Returns
        -------
        stacked : DataArray
            DataArray with stacked data.

        Examples
        --------
        >>> arr = xr.DataArray(
        ...     np.arange(6).reshape(2, 3),
        ...     coords=[("x", ["a", "b"]), ("y", [0, 1, 2])],
        ... )
        >>> arr
        <xarray.DataArray (x: 2, y: 3)> Size: 48B
        array([[0, 1, 2],
               [3, 4, 5]])
        Coordinates:
          * x        (x) <U1 8B 'a' 'b'
          * y        (y) int64 24B 0 1 2
        >>> stacked = arr.stack(z=("x", "y"))
        >>> stacked.indexes["z"]
        MultiIndex([('a', 0),
                    ('a', 1),
                    ('a', 2),
                    ('b', 0),
                    ('b', 1),
                    ('b', 2)],
                   name='z')

        See Also
        --------
        DataArray.unstack
        """
        ds = self._to_temp_dataset().stack(
            dim,
            create_index=create_index,
            index_cls=index_cls,
            **dim_kwargs,
        )
        return self._from_temp_dataset(ds)

    def unstack(
        self,
        dim: Dims = None,
        *,
        fill_value: Any = dtypes.NA,
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
            Value to be filled. If a dict-like, maps variable names to
            fill values. Use the data array's name to refer to its
            name. If not provided or if the dict-like does not contain
            all variables, the dtype's NA value will be used.
        sparse : bool, default: False
            Use sparse-array if True

        Returns
        -------
        unstacked : DataArray
            Array with unstacked data.

        Examples
        --------
        >>> arr = xr.DataArray(
        ...     np.arange(6).reshape(2, 3),
        ...     coords=[("x", ["a", "b"]), ("y", [0, 1, 2])],
        ... )
        >>> arr
        <xarray.DataArray (x: 2, y: 3)> Size: 48B
        array([[0, 1, 2],
               [3, 4, 5]])
        Coordinates:
          * x        (x) <U1 8B 'a' 'b'
          * y        (y) int64 24B 0 1 2
        >>> stacked = arr.stack(z=("x", "y"))
        >>> stacked.indexes["z"]
        MultiIndex([('a', 0),
                    ('a', 1),
                    ('a', 2),
                    ('b', 0),
                    ('b', 1),
                    ('b', 2)],
                   name='z')
        >>> roundtripped = stacked.unstack()
        >>> arr.identical(roundtripped)
        True

        See Also
        --------
        DataArray.stack
        """
        ds = self._to_temp_dataset().unstack(dim, fill_value=fill_value, sparse=sparse)
        return self._from_temp_dataset(ds)

    def to_unstacked_dataset(self, dim: Hashable, level: int | Hashable = 0) -> Dataset:
        """Unstack DataArray expanding to Dataset along a given level of a
        stacked coordinate.

        This is the inverse operation of Dataset.to_stacked_array.

        Parameters
        ----------
        dim : Hashable
            Name of existing dimension to unstack
        level : int or Hashable, default: 0
            The MultiIndex level to expand to a dataset along. Can either be
            the integer index of the level or its name.

        Returns
        -------
        unstacked: Dataset

        Examples
        --------
        >>> arr = xr.DataArray(
        ...     np.arange(6).reshape(2, 3),
        ...     coords=[("x", ["a", "b"]), ("y", [0, 1, 2])],
        ... )
        >>> data = xr.Dataset({"a": arr, "b": arr.isel(y=0)})
        >>> data
        <xarray.Dataset> Size: 96B
        Dimensions:  (x: 2, y: 3)
        Coordinates:
          * x        (x) <U1 8B 'a' 'b'
          * y        (y) int64 24B 0 1 2
        Data variables:
            a        (x, y) int64 48B 0 1 2 3 4 5
            b        (x) int64 16B 0 3
        >>> stacked = data.to_stacked_array("z", ["x"])
        >>> stacked.indexes["z"]
        MultiIndex([('a',   0),
                    ('a',   1),
                    ('a',   2),
                    ('b', nan)],
                   name='z')
        >>> roundtripped = stacked.to_unstacked_dataset(dim="z")
        >>> data.identical(roundtripped)
        True

        See Also
        --------
        Dataset.to_stacked_array
        """
        idx = self._indexes[dim].to_pandas_index()
        if not isinstance(idx, pd.MultiIndex):
            raise ValueError(f"'{dim}' is not a stacked coordinate")

        level_number = idx._get_level_number(level)  # type: ignore[attr-defined]
        variables = idx.levels[level_number]
        variable_dim = idx.names[level_number]

        # pull variables out of datarray
        data_dict = {}
        for k in variables:
            data_dict[k] = self.sel({variable_dim: k}, drop=True).squeeze(drop=True)

        # unstacked dataset
        return Dataset(data_dict)

    @deprecate_dims
    def transpose(
        self,
        *dim: Hashable,
        transpose_coords: bool = True,
        missing_dims: ErrorOptionsWithWarn = "raise",
    ) -> Self:
        """Return a new DataArray object with transposed dimensions.

        Parameters
        ----------
        *dim : Hashable, optional
            By default, reverse the dimensions. Otherwise, reorder the
            dimensions to this order.
        transpose_coords : bool, default: True
            If True, also transpose the coordinates of this DataArray.
        missing_dims : {"raise", "warn", "ignore"}, default: "raise"
            What to do if dimensions that should be selected from are not present in the
            DataArray:
            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        Returns
        -------
        transposed : DataArray
            The returned DataArray's array is transposed.

        Notes
        -----
        This operation returns a view of this array's data. It is
        lazy for dask-backed DataArrays but not for numpy-backed DataArrays
        -- the data will be fully loaded.

        See Also
        --------
        numpy.transpose
        Dataset.transpose
        """
        if dim:
            dim = tuple(infix_dims(dim, self.dims, missing_dims))
        variable = self.variable.transpose(*dim)
        if transpose_coords:
            coords: dict[Hashable, Variable] = {}
            for name, coord in self.coords.items():
                coord_dims = tuple(d for d in dim if d in coord.dims)
                coords[name] = coord.variable.transpose(*coord_dims)
            return self._replace(variable, coords)
        else:
            return self._replace(variable)

    @property
    def T(self) -> Self:
        return self.transpose()

    def drop_vars(
        self,
        names: str | Iterable[Hashable] | Callable[[Self], str | Iterable[Hashable]],
        *,
        errors: ErrorOptions = "raise",
    ) -> Self:
        """Returns an array with dropped variables.

        Parameters
        ----------
        names : Hashable or iterable of Hashable or Callable
            Name(s) of variables to drop. If a Callable, this object is passed as its
            only argument and its result is used.
        errors : {"raise", "ignore"}, default: "raise"
            If 'raise', raises a ValueError error if any of the variable
            passed are not in the dataset. If 'ignore', any given names that are in the
            DataArray are dropped and no error is raised.

        Returns
        -------
        dropped : Dataset
            New Dataset copied from `self` with variables removed.

        Examples
        -------
        >>> data = np.arange(12).reshape(4, 3)
        >>> da = xr.DataArray(
        ...     data=data,
        ...     dims=["x", "y"],
        ...     coords={"x": [10, 20, 30, 40], "y": [70, 80, 90]},
        ... )
        >>> da
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * x        (x) int64 32B 10 20 30 40
          * y        (y) int64 24B 70 80 90

        Removing a single variable:

        >>> da.drop_vars("x")
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * y        (y) int64 24B 70 80 90
        Dimensions without coordinates: x

        Removing a list of variables:

        >>> da.drop_vars(["x", "y"])
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Dimensions without coordinates: x, y

        >>> da.drop_vars(lambda x: x.coords)
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Dimensions without coordinates: x, y
        """
        if callable(names):
            names = names(self)
        ds = self._to_temp_dataset().drop_vars(names, errors=errors)
        return self._from_temp_dataset(ds)

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
        dropped : DataArray
            A new dataarray with dropped indexes.
        """
        ds = self._to_temp_dataset().drop_indexes(coord_names, errors=errors)
        return self._from_temp_dataset(ds)

    def drop(
        self,
        labels: Mapping[Any, Any] | None = None,
        dim: Hashable | None = None,
        *,
        errors: ErrorOptions = "raise",
        **labels_kwargs,
    ) -> Self:
        """Backward compatible method based on `drop_vars` and `drop_sel`

        Using either `drop_vars` or `drop_sel` is encouraged

        See Also
        --------
        DataArray.drop_vars
        DataArray.drop_sel
        """
        ds = self._to_temp_dataset().drop(labels, dim, errors=errors, **labels_kwargs)
        return self._from_temp_dataset(ds)

    def drop_sel(
        self,
        labels: Mapping[Any, Any] | None = None,
        *,
        errors: ErrorOptions = "raise",
        **labels_kwargs,
    ) -> Self:
        """Drop index labels from this DataArray.

        Parameters
        ----------
        labels : mapping of Hashable to Any
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
        dropped : DataArray

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.arange(25).reshape(5, 5),
        ...     coords={"x": np.arange(0, 9, 2), "y": np.arange(0, 13, 3)},
        ...     dims=("x", "y"),
        ... )
        >>> da
        <xarray.DataArray (x: 5, y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Coordinates:
          * x        (x) int64 40B 0 2 4 6 8
          * y        (y) int64 40B 0 3 6 9 12

        >>> da.drop_sel(x=[0, 2], y=9)
        <xarray.DataArray (x: 3, y: 4)> Size: 96B
        array([[10, 11, 12, 14],
               [15, 16, 17, 19],
               [20, 21, 22, 24]])
        Coordinates:
          * x        (x) int64 24B 4 6 8
          * y        (y) int64 32B 0 3 6 12

        >>> da.drop_sel({"x": 6, "y": [0, 3]})
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 2,  3,  4],
               [ 7,  8,  9],
               [12, 13, 14],
               [22, 23, 24]])
        Coordinates:
          * x        (x) int64 32B 0 2 4 8
          * y        (y) int64 24B 6 9 12
        """
        if labels_kwargs or isinstance(labels, dict):
            labels = either_dict_or_kwargs(labels, labels_kwargs, "drop")

        ds = self._to_temp_dataset().drop_sel(labels, errors=errors)
        return self._from_temp_dataset(ds)

    def drop_isel(
        self, indexers: Mapping[Any, Any] | None = None, **indexers_kwargs
    ) -> Self:
        """Drop index positions from this DataArray.

        Parameters
        ----------
        indexers : mapping of Hashable to Any or None, default: None
            Index locations to drop
        **indexers_kwargs : {dim: position, ...}, optional
            The keyword arguments form of ``dim`` and ``positions``

        Returns
        -------
        dropped : DataArray

        Raises
        ------
        IndexError

        Examples
        --------
        >>> da = xr.DataArray(np.arange(25).reshape(5, 5), dims=("X", "Y"))
        >>> da
        <xarray.DataArray (X: 5, Y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Dimensions without coordinates: X, Y

        >>> da.drop_isel(X=[0, 4], Y=2)
        <xarray.DataArray (X: 3, Y: 4)> Size: 96B
        array([[ 5,  6,  8,  9],
               [10, 11, 13, 14],
               [15, 16, 18, 19]])
        Dimensions without coordinates: X, Y

        >>> da.drop_isel({"X": 3, "Y": 3})
        <xarray.DataArray (X: 4, Y: 4)> Size: 128B
        array([[ 0,  1,  2,  4],
               [ 5,  6,  7,  9],
               [10, 11, 12, 14],
               [20, 21, 22, 24]])
        Dimensions without coordinates: X, Y
        """
        dataset = self._to_temp_dataset()
        dataset = dataset.drop_isel(indexers=indexers, **indexers_kwargs)
        return self._from_temp_dataset(dataset)

    def dropna(
        self,
        dim: Hashable,
        *,
        how: Literal["any", "all"] = "any",
        thresh: int | None = None,
    ) -> Self:
        """Returns a new array with dropped labels for missing values along
        the provided dimension.

        Parameters
        ----------
        dim : Hashable
            Dimension along which to drop missing values. Dropping along
            multiple dimensions simultaneously is not yet supported.
        how : {"any", "all"}, default: "any"
            - any : if any NA values are present, drop that label
            - all : if all values are NA, drop that label

        thresh : int or None, default: None
            If supplied, require this many non-NA values.

        Returns
        -------
        dropped : DataArray

        Examples
        --------
        >>> temperature = [
        ...     [0, 4, 2, 9],
        ...     [np.nan, np.nan, np.nan, np.nan],
        ...     [np.nan, 4, 2, 0],
        ...     [3, 1, 0, 0],
        ... ]
        >>> da = xr.DataArray(
        ...     data=temperature,
        ...     dims=["Y", "X"],
        ...     coords=dict(
        ...         lat=("Y", np.array([-20.0, -20.25, -20.50, -20.75])),
        ...         lon=("X", np.array([10.0, 10.25, 10.5, 10.75])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (Y: 4, X: 4)> Size: 128B
        array([[ 0.,  4.,  2.,  9.],
               [nan, nan, nan, nan],
               [nan,  4.,  2.,  0.],
               [ 3.,  1.,  0.,  0.]])
        Coordinates:
            lat      (Y) float64 32B -20.0 -20.25 -20.5 -20.75
            lon      (X) float64 32B 10.0 10.25 10.5 10.75
        Dimensions without coordinates: Y, X

        >>> da.dropna(dim="Y", how="any")
        <xarray.DataArray (Y: 2, X: 4)> Size: 64B
        array([[0., 4., 2., 9.],
               [3., 1., 0., 0.]])
        Coordinates:
            lat      (Y) float64 16B -20.0 -20.75
            lon      (X) float64 32B 10.0 10.25 10.5 10.75
        Dimensions without coordinates: Y, X

        Drop values only if all values along the dimension are NaN:

        >>> da.dropna(dim="Y", how="all")
        <xarray.DataArray (Y: 3, X: 4)> Size: 96B
        array([[ 0.,  4.,  2.,  9.],
               [nan,  4.,  2.,  0.],
               [ 3.,  1.,  0.,  0.]])
        Coordinates:
            lat      (Y) float64 24B -20.0 -20.5 -20.75
            lon      (X) float64 32B 10.0 10.25 10.5 10.75
        Dimensions without coordinates: Y, X
        """
        ds = self._to_temp_dataset().dropna(dim, how=how, thresh=thresh)
        return self._from_temp_dataset(ds)

    def fillna(self, value: Any) -> Self:
        """Fill missing values in this object.

        This operation follows the normal broadcasting and alignment rules that
        xarray uses for binary arithmetic, except the result is aligned to this
        object (``join='left'``) instead of aligned to the intersection of
        index coordinates (``join='inner'``).

        Parameters
        ----------
        value : scalar, ndarray or DataArray
            Used to fill all matching missing values in this array. If the
            argument is a DataArray, it is first aligned with (reindexed to)
            this array.

        Returns
        -------
        filled : DataArray

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 4, np.nan, 0, 3, np.nan]),
        ...     dims="Z",
        ...     coords=dict(
        ...         Z=("Z", np.arange(6)),
        ...         height=("Z", np.array([0, 10, 20, 30, 40, 50])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (Z: 6)> Size: 48B
        array([ 1.,  4., nan,  0.,  3., nan])
        Coordinates:
          * Z        (Z) int64 48B 0 1 2 3 4 5
            height   (Z) int64 48B 0 10 20 30 40 50

        Fill all NaN values with 0:

        >>> da.fillna(0)
        <xarray.DataArray (Z: 6)> Size: 48B
        array([1., 4., 0., 0., 3., 0.])
        Coordinates:
          * Z        (Z) int64 48B 0 1 2 3 4 5
            height   (Z) int64 48B 0 10 20 30 40 50

        Fill NaN values with corresponding values in array:

        >>> da.fillna(np.array([2, 9, 4, 2, 8, 9]))
        <xarray.DataArray (Z: 6)> Size: 48B
        array([1., 4., 4., 0., 3., 9.])
        Coordinates:
          * Z        (Z) int64 48B 0 1 2 3 4 5
            height   (Z) int64 48B 0 10 20 30 40 50
        """
        if utils.is_dict_like(value):
            raise TypeError(
                "cannot provide fill value as a dictionary with fillna on a DataArray"
            )
        out = ops.fillna(self, value)
        return out

    def interpolate_na(
        self,
        dim: Hashable | None = None,
        method: InterpOptions = "linear",
        limit: int | None = None,
        use_coordinate: bool | str = True,
        max_gap: (
            None
            | int
            | float
            | str
            | pd.Timedelta
            | np.timedelta64
            | datetime.timedelta
        ) = None,
        keep_attrs: bool | None = None,
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

        use_coordinate : bool or str, default: True
            Specifies which index to use as the x values in the interpolation
            formulated as `y = f(x)`. If False, values are treated as if
            equally-spaced along ``dim``. If True, the IndexVariable `dim` is
            used. If ``use_coordinate`` is a string, it specifies the name of a
            coordinate variable to use as the index.
        limit : int or None, default: None
            Maximum number of consecutive NaNs to fill. Must be greater than 0
            or None for no limit. This filling is done regardless of the size of
            the gap in the data. To only interpolate over gaps less than a given length,
            see ``max_gap``.
        max_gap : int, float, str, pandas.Timedelta, numpy.timedelta64, datetime.timedelta, default: None
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
        keep_attrs : bool or None, default: None
            If True, the dataarray's attributes (`attrs`) will be copied from
            the original object to the new one.  If False, the new
            object will be returned without attributes.
        **kwargs : dict, optional
            parameters passed verbatim to the underlying interpolation function

        Returns
        -------
        interpolated: DataArray
            Filled in DataArray.

        See Also
        --------
        numpy.interp
        scipy.interpolate

        Examples
        --------
        >>> da = xr.DataArray(
        ...     [np.nan, 2, 3, np.nan, 0], dims="x", coords={"x": [0, 1, 2, 3, 4]}
        ... )
        >>> da
        <xarray.DataArray (x: 5)> Size: 40B
        array([nan,  2.,  3., nan,  0.])
        Coordinates:
          * x        (x) int64 40B 0 1 2 3 4

        >>> da.interpolate_na(dim="x", method="linear")
        <xarray.DataArray (x: 5)> Size: 40B
        array([nan, 2. , 3. , 1.5, 0. ])
        Coordinates:
          * x        (x) int64 40B 0 1 2 3 4

        >>> da.interpolate_na(dim="x", method="linear", fill_value="extrapolate")
        <xarray.DataArray (x: 5)> Size: 40B
        array([1. , 2. , 3. , 1.5, 0. ])
        Coordinates:
          * x        (x) int64 40B 0 1 2 3 4
        """
        from xarray.core.missing import interp_na

        return interp_na(
            self,
            dim=dim,
            method=method,
            limit=limit,
            use_coordinate=use_coordinate,
            max_gap=max_gap,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def ffill(self, dim: Hashable, limit: int | None = None) -> Self:
        """Fill NaN values by propagating values forward

        *Requires bottleneck.*

        Parameters
        ----------
        dim : Hashable
            Specifies the dimension along which to propagate values when
            filling.
        limit : int or None, default: None
            The maximum number of consecutive NaN values to forward fill. In
            other words, if there is a gap with more than this number of
            consecutive NaNs, it will only be partially filled. Must be greater
            than 0 or None for no limit. Must be None or greater than or equal
            to axis length if filling along chunked axes (dimensions).

        Returns
        -------
        filled : DataArray

        Examples
        --------
        >>> temperature = np.array(
        ...     [
        ...         [np.nan, 1, 3],
        ...         [0, np.nan, 5],
        ...         [5, np.nan, np.nan],
        ...         [3, np.nan, np.nan],
        ...         [0, 2, 0],
        ...     ]
        ... )
        >>> da = xr.DataArray(
        ...     data=temperature,
        ...     dims=["Y", "X"],
        ...     coords=dict(
        ...         lat=("Y", np.array([-20.0, -20.25, -20.50, -20.75, -21.0])),
        ...         lon=("X", np.array([10.0, 10.25, 10.5])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (Y: 5, X: 3)> Size: 120B
        array([[nan,  1.,  3.],
               [ 0., nan,  5.],
               [ 5., nan, nan],
               [ 3., nan, nan],
               [ 0.,  2.,  0.]])
        Coordinates:
            lat      (Y) float64 40B -20.0 -20.25 -20.5 -20.75 -21.0
            lon      (X) float64 24B 10.0 10.25 10.5
        Dimensions without coordinates: Y, X

        Fill all NaN values:

        >>> da.ffill(dim="Y", limit=None)
        <xarray.DataArray (Y: 5, X: 3)> Size: 120B
        array([[nan,  1.,  3.],
               [ 0.,  1.,  5.],
               [ 5.,  1.,  5.],
               [ 3.,  1.,  5.],
               [ 0.,  2.,  0.]])
        Coordinates:
            lat      (Y) float64 40B -20.0 -20.25 -20.5 -20.75 -21.0
            lon      (X) float64 24B 10.0 10.25 10.5
        Dimensions without coordinates: Y, X

        Fill only the first of consecutive NaN values:

        >>> da.ffill(dim="Y", limit=1)
        <xarray.DataArray (Y: 5, X: 3)> Size: 120B
        array([[nan,  1.,  3.],
               [ 0.,  1.,  5.],
               [ 5., nan,  5.],
               [ 3., nan, nan],
               [ 0.,  2.,  0.]])
        Coordinates:
            lat      (Y) float64 40B -20.0 -20.25 -20.5 -20.75 -21.0
            lon      (X) float64 24B 10.0 10.25 10.5
        Dimensions without coordinates: Y, X
        """
        from xarray.core.missing import ffill

        return ffill(self, dim, limit=limit)

    def bfill(self, dim: Hashable, limit: int | None = None) -> Self:
        """Fill NaN values by propagating values backward

        *Requires bottleneck.*

        Parameters
        ----------
        dim : str
            Specifies the dimension along which to propagate values when
            filling.
        limit : int or None, default: None
            The maximum number of consecutive NaN values to backward fill. In
            other words, if there is a gap with more than this number of
            consecutive NaNs, it will only be partially filled. Must be greater
            than 0 or None for no limit. Must be None or greater than or equal
            to axis length if filling along chunked axes (dimensions).

        Returns
        -------
        filled : DataArray

        Examples
        --------
        >>> temperature = np.array(
        ...     [
        ...         [0, 1, 3],
        ...         [0, np.nan, 5],
        ...         [5, np.nan, np.nan],
        ...         [3, np.nan, np.nan],
        ...         [np.nan, 2, 0],
        ...     ]
        ... )
        >>> da = xr.DataArray(
        ...     data=temperature,
        ...     dims=["Y", "X"],
        ...     coords=dict(
        ...         lat=("Y", np.array([-20.0, -20.25, -20.50, -20.75, -21.0])),
        ...         lon=("X", np.array([10.0, 10.25, 10.5])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (Y: 5, X: 3)> Size: 120B
        array([[ 0.,  1.,  3.],
               [ 0., nan,  5.],
               [ 5., nan, nan],
               [ 3., nan, nan],
               [nan,  2.,  0.]])
        Coordinates:
            lat      (Y) float64 40B -20.0 -20.25 -20.5 -20.75 -21.0
            lon      (X) float64 24B 10.0 10.25 10.5
        Dimensions without coordinates: Y, X

        Fill all NaN values:

        >>> da.bfill(dim="Y", limit=None)
        <xarray.DataArray (Y: 5, X: 3)> Size: 120B
        array([[ 0.,  1.,  3.],
               [ 0.,  2.,  5.],
               [ 5.,  2.,  0.],
               [ 3.,  2.,  0.],
               [nan,  2.,  0.]])
        Coordinates:
            lat      (Y) float64 40B -20.0 -20.25 -20.5 -20.75 -21.0
            lon      (X) float64 24B 10.0 10.25 10.5
        Dimensions without coordinates: Y, X

        Fill only the first of consecutive NaN values:

        >>> da.bfill(dim="Y", limit=1)
        <xarray.DataArray (Y: 5, X: 3)> Size: 120B
        array([[ 0.,  1.,  3.],
               [ 0., nan,  5.],
               [ 5., nan, nan],
               [ 3.,  2.,  0.],
               [nan,  2.,  0.]])
        Coordinates:
            lat      (Y) float64 40B -20.0 -20.25 -20.5 -20.75 -21.0
            lon      (X) float64 24B 10.0 10.25 10.5
        Dimensions without coordinates: Y, X
        """
        from xarray.core.missing import bfill

        return bfill(self, dim, limit=limit)

    def combine_first(self, other: Self) -> Self:
        """Combine two DataArray objects, with union of coordinates.

        This operation follows the normal broadcasting and alignment rules of
        ``join='outer'``.  Default to non-null values of array calling the
        method.  Use np.nan to fill in vacant cells after alignment.

        Parameters
        ----------
        other : DataArray
            Used to fill all matching missing values in this array.

        Returns
        -------
        DataArray
        """
        return ops.fillna(self, other, join="outer")

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Self:
        """Reduce this array by applying `func` along some dimension(s).

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `f(x, axis=axis, **kwargs)` to return the result of reducing an
            np.ndarray over an integer valued axis.
        dim : "...", str, Iterable of Hashable or None, optional
            Dimension(s) over which to apply `func`. By default `func` is
            applied over all dimensions.
        axis : int or sequence of int, optional
            Axis(es) over which to repeatedly apply `func`. Only one of the
            'dim' and 'axis' arguments can be supplied. If neither are
            supplied, then the reduction is calculated over the flattened array
            (by calling `f(x)` without an axis argument).
        keep_attrs : bool or None, optional
            If True, the variable's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        keepdims : bool, default: False
            If True, the dimensions which are reduced are left in the result
            as dimensions of size one. Coordinates that use these dimensions
            are removed.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : DataArray
            DataArray with this object's array replaced with an array with
            summarized data and the indicated dimension(s) removed.
        """

        var = self.variable.reduce(func, dim, axis, keep_attrs, keepdims, **kwargs)
        return self._replace_maybe_drop_dims(var)

    def to_pandas(self) -> Self | pd.Series | pd.DataFrame:
        """Convert this array into a pandas object with the same shape.

        The type of the returned object depends on the number of DataArray
        dimensions:

        * 0D -> `xarray.DataArray`
        * 1D -> `pandas.Series`
        * 2D -> `pandas.DataFrame`

        Only works for arrays with 2 or fewer dimensions.

        The DataArray constructor performs the inverse transformation.

        Returns
        -------
        result : DataArray | Series | DataFrame
            DataArray, pandas Series or pandas DataFrame.
        """
        # TODO: consolidate the info about pandas constructors and the
        # attributes that correspond to their indexes into a separate module?
        constructors: dict[int, Callable] = {
            0: lambda x: x,
            1: pd.Series,
            2: pd.DataFrame,
        }
        try:
            constructor = constructors[self.ndim]
        except KeyError as err:
            raise ValueError(
                f"Cannot convert arrays with {self.ndim} dimensions into "
                "pandas objects. Requires 2 or fewer dimensions."
            ) from err
        indexes = [self.get_index(dim) for dim in self.dims]
        if isinstance(self._variable._data, PandasExtensionArray):
            values = self._variable._data.array
        else:
            values = self.values
        pandas_object = constructor(values, *indexes)
        if isinstance(pandas_object, pd.Series):
            pandas_object.name = self.name
        return pandas_object

    def to_dataframe(
        self, name: Hashable | None = None, dim_order: Sequence[Hashable] | None = None
    ) -> pd.DataFrame:
        """Convert this array and its coordinates into a tidy pandas.DataFrame.

        The DataFrame is indexed by the Cartesian product of index coordinates
        (in the form of a :py:class:`pandas.MultiIndex`). Other coordinates are
        included as columns in the DataFrame.

        For 1D and 2D DataArrays, see also :py:func:`DataArray.to_pandas` which
        doesn't rely on a MultiIndex to build the DataFrame.

        Parameters
        ----------
        name: Hashable or None, optional
            Name to give to this array (required if unnamed).
        dim_order: Sequence of Hashable or None, optional
            Hierarchical dimension order for the resulting dataframe.
            Array content is transposed to this order and then written out as flat
            vectors in contiguous order, so the last dimension in this list
            will be contiguous in the resulting DataFrame. This has a major
            influence on which operations are efficient on the resulting
            dataframe.

            If provided, must include all dimensions of this DataArray. By default,
            dimensions are sorted according to the DataArray dimensions order.

        Returns
        -------
        result: DataFrame
            DataArray as a pandas DataFrame.

        See also
        --------
        DataArray.to_pandas
        DataArray.to_series
        """
        if name is None:
            name = self.name
        if name is None:
            raise ValueError(
                "cannot convert an unnamed DataArray to a "
                "DataFrame: use the ``name`` parameter"
            )
        if self.ndim == 0:
            raise ValueError("cannot convert a scalar to a DataFrame")

        # By using a unique name, we can convert a DataArray into a DataFrame
        # even if it shares a name with one of its coordinates.
        # I would normally use unique_name = object() but that results in a
        # dataframe with columns in the wrong order, for reasons I have not
        # been able to debug (possibly a pandas bug?).
        unique_name = "__unique_name_identifier_z98xfz98xugfg73ho__"
        ds = self._to_dataset_whole(name=unique_name)

        if dim_order is None:
            ordered_dims = dict(zip(self.dims, self.shape, strict=True))
        else:
            ordered_dims = ds._normalize_dim_order(dim_order=dim_order)

        df = ds._to_dataframe(ordered_dims)
        df.columns = [name if c == unique_name else c for c in df.columns]
        return df

    def to_series(self) -> pd.Series:
        """Convert this array into a pandas.Series.

        The Series is indexed by the Cartesian product of index coordinates
        (in the form of a :py:class:`pandas.MultiIndex`).

        Returns
        -------
        result : Series
            DataArray as a pandas Series.

        See also
        --------
        DataArray.to_pandas
        DataArray.to_dataframe
        """
        index = self.coords.to_index()
        return pd.Series(self.values.reshape(-1), index=index, name=self.name)

    def to_masked_array(self, copy: bool = True) -> np.ma.MaskedArray:
        """Convert this array into a numpy.ma.MaskedArray

        Parameters
        ----------
        copy : bool, default: True
            If True make a copy of the array in the result. If False,
            a MaskedArray view of DataArray.values is returned.

        Returns
        -------
        result : MaskedArray
            Masked where invalid values (nan or inf) occur.
        """
        values = self.to_numpy()  # only compute lazy arrays once
        isnull = pd.isnull(values)
        return np.ma.MaskedArray(data=values, mask=isnull, copy=copy)

    # path=None writes to bytes
    @overload
    def to_netcdf(
        self,
        path: None = None,
        mode: NetcdfWriteModes = "w",
        format: T_NetcdfTypes | None = None,
        group: str | None = None,
        engine: T_NetcdfEngine | None = None,
        encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
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
        encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
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
        path: str | PathLike,
        mode: NetcdfWriteModes = "w",
        format: T_NetcdfTypes | None = None,
        group: str | None = None,
        engine: T_NetcdfEngine | None = None,
        encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
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
        encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
        unlimited_dims: Iterable[Hashable] | None = None,
        compute: bool = True,
        invalid_netcdf: bool = False,
        auto_complex: bool | None = None,
    ) -> Delayed | None: ...

    def to_netcdf(
        self,
        path: str | PathLike | None = None,
        mode: NetcdfWriteModes = "w",
        format: T_NetcdfTypes | None = None,
        group: str | None = None,
        engine: T_NetcdfEngine | None = None,
        encoding: Mapping[Hashable, Mapping[str, Any]] | None = None,
        unlimited_dims: Iterable[Hashable] | None = None,
        compute: bool = True,
        invalid_netcdf: bool = False,
        auto_complex: bool | None = None,
    ) -> memoryview | Delayed | None:
        """Write DataArray contents to a netCDF file.

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
            "zlib": True}, ...}``

            The `h5netcdf` engine supports both the NetCDF4-style compression
            encoding parameters ``{"zlib": True, "complevel": 9}`` and the h5py
            ones ``{"compression": "gzip", "compression_opts": 9}``.
            This allows using any compression plugin installed in the HDF5
            library, e.g. LZF.

        unlimited_dims : iterable of Hashable, optional
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
            * None otherwise

        Notes
        -----
        Only xarray.Dataset objects can be written to netCDF files, so
        the xarray.DataArray is converted to a xarray.Dataset object
        containing a single variable. If the DataArray has no name, or if the
        name is the same as a coordinate name, then it is given the name
        ``"__xarray_dataarray_variable__"``.

        [netCDF4 backend only] netCDF4 enums are decoded into the
        dataarray dtype metadata.

        See Also
        --------
        Dataset.to_netcdf
        """
        from xarray.backends.api import DATAARRAY_NAME, DATAARRAY_VARIABLE, to_netcdf

        if self.name is None:
            # If no name is set then use a generic xarray name
            dataset = self.to_dataset(name=DATAARRAY_VARIABLE)
        elif self.name in self.coords or self.name in self.dims:
            # The name is the same as one of the coords names, which netCDF
            # doesn't support, so rename it but keep track of the old name
            dataset = self.to_dataset(name=DATAARRAY_VARIABLE)
            dataset.attrs[DATAARRAY_NAME] = self.name
        else:
            # No problems with the name - so we're fine!
            dataset = self.to_dataset()

        return to_netcdf(  # type: ignore[return-value]  # mypy cannot resolve the overloads:(
            dataset,
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
        *,
        encoding: Mapping | None = None,
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
        """Write DataArray contents to a Zarr store

        Zarr chunks are determined in the following way:

        - From the ``chunks`` attribute in each variable's ``encoding``
          (can be set via `DataArray.chunk`).
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
            If True, apply zarr's `consolidate_metadata` function to the store
            after writing metadata and read existing stores with consolidated
            metadata; if False, do not. The default (`consolidated=None`) means
            write consolidated metadata and attempt to read consolidated
            metadata for existing stores (falling back to non-consolidated).

            When the experimental ``zarr_version=3``, ``consolidated`` must be
            either be ``None`` or ``False``.
        append_dim : hashable, optional
            If set, the dimension along which the data will be appended. All
            other dimensions on overridden variables must remain the same size.
        region : dict, optional
            Optional mapping from dimension names to integer slices along
            dataarray dimensions to indicate the region of existing zarr array(s)
            in which to write this datarray's data. For example,
            ``{'x': slice(0, 1000), 'y': slice(10000, 11000)}`` would indicate
            that values should be written to the region ``0:1000`` along ``x``
            and ``10000:11000`` along ``y``.

            Two restrictions apply to the use of ``region``:

            - If ``region`` is set, _all_ variables in a dataarray must have at
              least one dimension in common with the region. Other variables
              should be written in a separate call to ``to_zarr()``.
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
            with ``compute=False`` to initialize a Zarr store from an existing
            DataArray with arbitrary chunk structure.
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
        Dataset.to_zarr
        :ref:`io.zarr`
            The I/O user guide, with more details and examples.
        """
        from xarray.backends.api import DATAARRAY_NAME, DATAARRAY_VARIABLE, to_zarr

        if self.name is None:
            # If no name is set then use a generic xarray name
            dataset = self.to_dataset(name=DATAARRAY_VARIABLE)
        elif self.name in self.coords or self.name in self.dims:
            # The name is the same as one of the coords names, which the netCDF data model
            # does not support, so rename it but keep track of the old name
            dataset = self.to_dataset(name=DATAARRAY_VARIABLE)
            dataset.attrs[DATAARRAY_NAME] = self.name
        else:
            # No problems with the name - so we're fine!
            dataset = self.to_dataset()

        return to_zarr(  # type: ignore[call-overload,misc]
            dataset,
            store=store,
            chunk_store=chunk_store,
            mode=mode,
            synchronizer=synchronizer,
            group=group,
            encoding=encoding,
            compute=compute,
            consolidated=consolidated,
            append_dim=append_dim,
            region=region,
            safe_chunks=safe_chunks,
            align_chunks=align_chunks,
            storage_options=storage_options,
            zarr_version=zarr_version,
            zarr_format=zarr_format,
            write_empty_chunks=write_empty_chunks,
            chunkmanager_store_kwargs=chunkmanager_store_kwargs,
        )

    def to_dict(
        self, data: bool | Literal["list", "array"] = "list", encoding: bool = False
    ) -> dict[str, Any]:
        """
        Convert this xarray.DataArray into a dictionary following xarray
        naming conventions.

        Converts all variables and attributes to native Python objects.
        Useful for converting to json. To avoid datetime incompatibility
        use decode_times=False kwarg in xarray.open_dataset.

        Parameters
        ----------
        data : bool or {"list", "array"}, default: "list"
            Whether to include the actual data in the dictionary. When set to
            False, returns just the schema. If set to "array", returns data as
            underlying array type. If set to "list" (or True for backwards
            compatibility), returns data in lists of Python data types. Note
            that for obtaining the "list" output efficiently, use
            `da.compute().to_dict(data="list")`.

        encoding : bool, default: False
            Whether to include the Dataset's encoding in the dictionary.

        Returns
        -------
        dict: dict

        See Also
        --------
        DataArray.from_dict
        Dataset.to_dict
        """
        d = self.variable.to_dict(data=data)
        d.update({"coords": {}, "name": self.name})
        for k, coord in self.coords.items():
            d["coords"][k] = coord.variable.to_dict(data=data)
        if encoding:
            d["encoding"] = dict(self.encoding)
        return d

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> Self:
        """Convert a dictionary into an xarray.DataArray

        Parameters
        ----------
        d : dict
            Mapping with a minimum structure of {"dims": [...], "data": [...]}

        Returns
        -------
        obj : xarray.DataArray

        See Also
        --------
        DataArray.to_dict
        Dataset.from_dict

        Examples
        --------
        >>> d = {"dims": "t", "data": [1, 2, 3]}
        >>> da = xr.DataArray.from_dict(d)
        >>> da
        <xarray.DataArray (t: 3)> Size: 24B
        array([1, 2, 3])
        Dimensions without coordinates: t

        >>> d = {
        ...     "coords": {
        ...         "t": {"dims": "t", "data": [0, 1, 2], "attrs": {"units": "s"}}
        ...     },
        ...     "attrs": {"title": "air temperature"},
        ...     "dims": "t",
        ...     "data": [10, 20, 30],
        ...     "name": "a",
        ... }
        >>> da = xr.DataArray.from_dict(d)
        >>> da
        <xarray.DataArray 'a' (t: 3)> Size: 24B
        array([10, 20, 30])
        Coordinates:
          * t        (t) int64 24B 0 1 2
        Attributes:
            title:    air temperature
        """
        coords = None
        if "coords" in d:
            try:
                coords = {
                    k: (v["dims"], v["data"], v.get("attrs"))
                    for k, v in d["coords"].items()
                }
            except KeyError as e:
                raise ValueError(
                    f"cannot convert dict when coords are missing the key '{e.args[0]}'"
                ) from e
        try:
            data = d["data"]
        except KeyError as err:
            raise ValueError("cannot convert dict without the key 'data''") from err
        else:
            obj = cls(data, coords, d.get("dims"), d.get("name"), d.get("attrs"))

        obj.encoding.update(d.get("encoding", {}))

        return obj

    @classmethod
    def from_series(cls, series: pd.Series, sparse: bool = False) -> DataArray:
        """Convert a pandas.Series into an xarray.DataArray.

        If the series's index is a MultiIndex, it will be expanded into a
        tensor product of one-dimensional coordinates (filling in missing
        values with NaN). Thus this operation should be the inverse of the
        `to_series` method.

        Parameters
        ----------
        series : Series
            Pandas Series object to convert.
        sparse : bool, default: False
            If sparse=True, creates a sparse array instead of a dense NumPy array.
            Requires the pydata/sparse package.

        See Also
        --------
        DataArray.to_series
        Dataset.from_dataframe
        """
        temp_name = "__temporary_name"
        df = pd.DataFrame({temp_name: series})
        ds = Dataset.from_dataframe(df, sparse=sparse)
        result = ds[temp_name]
        result.name = series.name
        return result

    def to_iris(self) -> iris_Cube:
        """Convert this array into a iris.cube.Cube"""
        from xarray.convert import to_iris

        return to_iris(self)

    @classmethod
    def from_iris(cls, cube: iris_Cube) -> Self:
        """Convert a iris.cube.Cube into an xarray.DataArray"""
        from xarray.convert import from_iris

        return from_iris(cube)

    def _all_compat(self, other: Self, compat_str: str) -> bool:
        """Helper function for equals, broadcast_equals, and identical"""

        def compat(x, y):
            return getattr(x.variable, compat_str)(y.variable)

        return utils.dict_equiv(self.coords, other.coords, compat=compat) and compat(
            self, other
        )

    def broadcast_equals(self, other: Self) -> bool:
        """Two DataArrays are broadcast equal if they are equal after
        broadcasting them against each other such that they have the same
        dimensions.

        Parameters
        ----------
        other : DataArray
            DataArray to compare to.

        Returns
        ----------
        equal : bool
            True if the two DataArrays are broadcast equal.

        See Also
        --------
        DataArray.equals
        DataArray.identical

        Examples
        --------
        >>> a = xr.DataArray([1, 2], dims="X")
        >>> b = xr.DataArray([[1, 1], [2, 2]], dims=["X", "Y"])
        >>> a
        <xarray.DataArray (X: 2)> Size: 16B
        array([1, 2])
        Dimensions without coordinates: X
        >>> b
        <xarray.DataArray (X: 2, Y: 2)> Size: 32B
        array([[1, 1],
               [2, 2]])
        Dimensions without coordinates: X, Y

        .equals returns True if two DataArrays have the same values, dimensions, and coordinates. .broadcast_equals returns True if the results of broadcasting two DataArrays against each other have the same values, dimensions, and coordinates.

        >>> a.equals(b)
        False
        >>> a2, b2 = xr.broadcast(a, b)
        >>> a2.equals(b2)
        True
        >>> a.broadcast_equals(b)
        True
        """
        try:
            return self._all_compat(other, "broadcast_equals")
        except (TypeError, AttributeError):
            return False

    def equals(self, other: Self) -> bool:
        """True if two DataArrays have the same dimensions, coordinates and
        values; otherwise False.

        DataArrays can still be equal (like pandas objects) if they have NaN
        values in the same locations.

        This method is necessary because `v1 == v2` for ``DataArray``
        does element-wise comparisons (like numpy.ndarrays).

        Parameters
        ----------
        other : DataArray
            DataArray to compare to.

        Returns
        ----------
        equal : bool
            True if the two DataArrays are equal.

        See Also
        --------
        DataArray.broadcast_equals
        DataArray.identical

        Examples
        --------
        >>> a = xr.DataArray([1, 2, 3], dims="X")
        >>> b = xr.DataArray([1, 2, 3], dims="X", attrs=dict(units="m"))
        >>> c = xr.DataArray([1, 2, 3], dims="Y")
        >>> d = xr.DataArray([3, 2, 1], dims="X")
        >>> a
        <xarray.DataArray (X: 3)> Size: 24B
        array([1, 2, 3])
        Dimensions without coordinates: X
        >>> b
        <xarray.DataArray (X: 3)> Size: 24B
        array([1, 2, 3])
        Dimensions without coordinates: X
        Attributes:
            units:    m
        >>> c
        <xarray.DataArray (Y: 3)> Size: 24B
        array([1, 2, 3])
        Dimensions without coordinates: Y
        >>> d
        <xarray.DataArray (X: 3)> Size: 24B
        array([3, 2, 1])
        Dimensions without coordinates: X

        >>> a.equals(b)
        True
        >>> a.equals(c)
        False
        >>> a.equals(d)
        False
        """
        try:
            return self._all_compat(other, "equals")
        except (TypeError, AttributeError):
            return False

    def identical(self, other: Self) -> bool:
        """Like equals, but also checks the array name and attributes, and
        attributes on all coordinates.

        Parameters
        ----------
        other : DataArray
            DataArray to compare to.

        Returns
        ----------
        equal : bool
            True if the two DataArrays are identical.

        See Also
        --------
        DataArray.broadcast_equals
        DataArray.equals

        Examples
        --------
        >>> a = xr.DataArray([1, 2, 3], dims="X", attrs=dict(units="m"), name="Width")
        >>> b = xr.DataArray([1, 2, 3], dims="X", attrs=dict(units="m"), name="Width")
        >>> c = xr.DataArray([1, 2, 3], dims="X", attrs=dict(units="ft"), name="Width")
        >>> a
        <xarray.DataArray 'Width' (X: 3)> Size: 24B
        array([1, 2, 3])
        Dimensions without coordinates: X
        Attributes:
            units:    m
        >>> b
        <xarray.DataArray 'Width' (X: 3)> Size: 24B
        array([1, 2, 3])
        Dimensions without coordinates: X
        Attributes:
            units:    m
        >>> c
        <xarray.DataArray 'Width' (X: 3)> Size: 24B
        array([1, 2, 3])
        Dimensions without coordinates: X
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
        """
        try:
            return self.name == other.name and self._all_compat(other, "identical")
        except (TypeError, AttributeError):
            return False

    def __array_wrap__(self, obj, context=None, return_scalar=False) -> Self:
        new_var = self.variable.__array_wrap__(obj, context, return_scalar)
        return self._replace(new_var)

    def __matmul__(self, obj: T_Xarray) -> T_Xarray:
        return self.dot(obj)

    def __rmatmul__(self, other: T_Xarray) -> T_Xarray:
        # currently somewhat duplicative, as only other DataArrays are
        # compatible with matmul
        return computation.dot(other, self)

    def _unary_op(self, f: Callable, *args, **kwargs) -> Self:
        keep_attrs = kwargs.pop("keep_attrs", None)
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")
            warnings.filterwarnings(
                "ignore", r"Mean of empty slice", category=RuntimeWarning
            )
            with np.errstate(all="ignore"):
                da = self.__array_wrap__(f(self.variable.data, *args, **kwargs))
            if keep_attrs:
                da.attrs = self.attrs
            return da

    def _binary_op(
        self, other: DaCompatible, f: Callable, reflexive: bool = False
    ) -> Self:
        from xarray.core.datatree import DataTree
        from xarray.core.groupby import GroupBy

        if isinstance(other, DataTree | Dataset | GroupBy):
            return NotImplemented
        if isinstance(other, DataArray):
            align_type = OPTIONS["arithmetic_join"]
            self, other = align(self, other, join=align_type, copy=False)
        other_variable_or_arraylike: DaCompatible = getattr(other, "variable", other)
        other_coords = getattr(other, "coords", None)

        variable = (
            f(self.variable, other_variable_or_arraylike)
            if not reflexive
            else f(other_variable_or_arraylike, self.variable)
        )
        coords, indexes = self.coords._merge_raw(other_coords, reflexive)
        name = result_name([self, other])

        return self._replace(variable, coords, name, indexes=indexes)

    def _inplace_binary_op(self, other: DaCompatible, f: Callable) -> Self:
        from xarray.core.groupby import GroupBy

        if isinstance(other, GroupBy):
            raise TypeError(
                "in-place operations between a DataArray and "
                "a grouped object are not permitted"
            )
        # n.b. we can't align other to self (with other.reindex_like(self))
        # because `other` may be converted into floats, which would cause
        # in-place arithmetic to fail unpredictably. Instead, we simply
        # don't support automatic alignment with in-place arithmetic.
        other_coords = getattr(other, "coords", None)
        other_variable = getattr(other, "variable", other)
        try:
            with self.coords._merge_inplace(other_coords):
                f(self.variable, other_variable)
        except MergeError as exc:
            raise MergeError(
                "Automatic alignment is not supported for in-place operations.\n"
                "Consider aligning the indices manually or using a not-in-place operation.\n"
                "See https://github.com/pydata/xarray/issues/3910 for more explanations."
            ) from exc
        return self

    def _copy_attrs_from(self, other: DataArray | Dataset | Variable) -> None:
        self.attrs = other.attrs

    plot = utils.UncachedAccessor(DataArrayPlotAccessor)

    def _title_for_slice(self, truncate: int = 50) -> str:
        """
        If the dataarray has 1 dimensional coordinates or comes from a slice
        we can show that info in the title

        Parameters
        ----------
        truncate : int, default: 50
            maximum number of characters for title

        Returns
        -------
        title : string
            Can be used for plot titles

        """
        one_dims = []
        for dim, coord in self.coords.items():
            if coord.size == 1:
                one_dims.append(
                    f"{dim} = {format_item(coord.values)}{_get_units_from_attrs(coord)}"
                )

        title = ", ".join(one_dims)
        if len(title) > truncate:
            title = title[: (truncate - 3)] + "..."

        return title

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
        difference : DataArray
            The n-th order finite difference of this object.

        Notes
        -----
        `n` matches numpy's behavior and is different from pandas' first argument named
        `periods`.

        Examples
        --------
        >>> arr = xr.DataArray([5, 5, 6, 6], [[1, 2, 3, 4]], ["x"])
        >>> arr.diff("x")
        <xarray.DataArray (x: 3)> Size: 24B
        array([0, 1, 0])
        Coordinates:
          * x        (x) int64 24B 2 3 4
        >>> arr.diff("x", 2)
        <xarray.DataArray (x: 2)> Size: 16B
        array([ 1, -1])
        Coordinates:
          * x        (x) int64 16B 3 4

        See Also
        --------
        DataArray.differentiate
        """
        ds = self._to_temp_dataset().diff(n=n, dim=dim, label=label)
        return self._from_temp_dataset(ds)

    def shift(
        self,
        shifts: Mapping[Any, int] | None = None,
        fill_value: Any = dtypes.NA,
        **shifts_kwargs: int,
    ) -> Self:
        """Shift this DataArray by an offset along one or more dimensions.

        Only the data is moved; coordinates stay in place. This is consistent
        with the behavior of ``shift`` in pandas.

        Values shifted from beyond array bounds will appear at one end of
        each dimension, which are filled according to `fill_value`. For periodic
        offsets instead see `roll`.

        Parameters
        ----------
        shifts : mapping of Hashable to int or None, optional
            Integer offset to shift along each of the given dimensions.
            Positive offsets shift to the right; negative offsets shift to the
            left.
        fill_value : scalar, optional
            Value to use for newly missing values
        **shifts_kwargs
            The keyword arguments form of ``shifts``.
            One of shifts or shifts_kwargs must be provided.

        Returns
        -------
        shifted : DataArray
            DataArray with the same coordinates and attributes but shifted
            data.

        See Also
        --------
        roll

        Examples
        --------
        >>> arr = xr.DataArray([5, 6, 7], dims="x")
        >>> arr.shift(x=1)
        <xarray.DataArray (x: 3)> Size: 24B
        array([nan,  5.,  6.])
        Dimensions without coordinates: x
        """
        variable = self.variable.shift(
            shifts=shifts, fill_value=fill_value, **shifts_kwargs
        )
        return self._replace(variable=variable)

    def roll(
        self,
        shifts: Mapping[Hashable, int] | None = None,
        roll_coords: bool = False,
        **shifts_kwargs: int,
    ) -> Self:
        """Roll this array by an offset along one or more dimensions.

        Unlike shift, roll treats the given dimensions as periodic, so will not
        create any missing values to be filled.

        Unlike shift, roll may rotate all variables, including coordinates
        if specified. The direction of rotation is consistent with
        :py:func:`numpy.roll`.

        Parameters
        ----------
        shifts : mapping of Hashable to int, optional
            Integer offset to rotate each of the given dimensions.
            Positive offsets roll to the right; negative offsets roll to the
            left.
        roll_coords : bool, default: False
            Indicates whether to roll the coordinates by the offset too.
        **shifts_kwargs : {dim: offset, ...}, optional
            The keyword arguments form of ``shifts``.
            One of shifts or shifts_kwargs must be provided.

        Returns
        -------
        rolled : DataArray
            DataArray with the same attributes but rolled data and coordinates.

        See Also
        --------
        shift

        Examples
        --------
        >>> arr = xr.DataArray([5, 6, 7], dims="x")
        >>> arr.roll(x=1)
        <xarray.DataArray (x: 3)> Size: 24B
        array([7, 5, 6])
        Dimensions without coordinates: x
        """
        ds = self._to_temp_dataset().roll(
            shifts=shifts, roll_coords=roll_coords, **shifts_kwargs
        )
        return self._from_temp_dataset(ds)

    @property
    def real(self) -> Self:
        """
        The real part of the array.

        See Also
        --------
        numpy.ndarray.real
        """
        return self._replace(self.variable.real)

    @property
    def imag(self) -> Self:
        """
        The imaginary part of the array.

        See Also
        --------
        numpy.ndarray.imag
        """
        return self._replace(self.variable.imag)

    @deprecate_dims
    def dot(
        self,
        other: T_Xarray,
        dim: Dims = None,
    ) -> T_Xarray:
        """Perform dot product of two DataArrays along their shared dims.

        Equivalent to taking taking tensordot over all shared dims.

        Parameters
        ----------
        other : DataArray
            The other array with which the dot product is performed.
        dim : ..., str, Iterable of Hashable or None, optional
            Which dimensions to sum over. Ellipsis (`...`) sums over all dimensions.
            If not specified, then all the common dimensions are summed over.

        Returns
        -------
        result : DataArray
            Array resulting from the dot product over all shared dimensions.

        See Also
        --------
        dot
        numpy.tensordot

        Examples
        --------
        >>> da_vals = np.arange(6 * 5 * 4).reshape((6, 5, 4))
        >>> da = xr.DataArray(da_vals, dims=["x", "y", "z"])
        >>> dm_vals = np.arange(4)
        >>> dm = xr.DataArray(dm_vals, dims=["z"])

        >>> dm.dims
        ('z',)

        >>> da.dims
        ('x', 'y', 'z')

        >>> dot_result = da.dot(dm)
        >>> dot_result.dims
        ('x', 'y')

        """
        if isinstance(other, Dataset):
            raise NotImplementedError(
                "dot products are not yet supported with Dataset objects."
            )
        if not isinstance(other, DataArray):
            raise TypeError("dot only operates on DataArrays.")

        return computation.dot(self, other, dim=dim)

    def sortby(
        self,
        variables: (
            Hashable
            | DataArray
            | Sequence[Hashable | DataArray]
            | Callable[[Self], Hashable | DataArray | Sequence[Hashable | DataArray]]
        ),
        ascending: bool = True,
    ) -> Self:
        """Sort object by labels or values (along an axis).

        Sorts the dataarray, either along specified dimensions,
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
        sorted : DataArray
            A new dataarray where all the specified dims are sorted by dim
            labels.

        See Also
        --------
        Dataset.sortby
        numpy.sort
        pandas.sort_values
        pandas.sort_index

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.arange(5, 0, -1),
        ...     coords=[pd.date_range("1/1/2000", periods=5)],
        ...     dims="time",
        ... )
        >>> da
        <xarray.DataArray (time: 5)> Size: 40B
        array([5, 4, 3, 2, 1])
        Coordinates:
          * time     (time) datetime64[ns] 40B 2000-01-01 2000-01-02 ... 2000-01-05

        >>> da.sortby(da)
        <xarray.DataArray (time: 5)> Size: 40B
        array([1, 2, 3, 4, 5])
        Coordinates:
          * time     (time) datetime64[ns] 40B 2000-01-05 2000-01-04 ... 2000-01-01

        >>> da.sortby(lambda x: x)
        <xarray.DataArray (time: 5)> Size: 40B
        array([1, 2, 3, 4, 5])
        Coordinates:
          * time     (time) datetime64[ns] 40B 2000-01-05 2000-01-04 ... 2000-01-01
        """
        # We need to convert the callable here rather than pass it through to the
        # dataset method, since otherwise the dataset method would try to call the
        # callable with the dataset as the object
        if callable(variables):
            variables = variables(self)
        ds = self._to_temp_dataset().sortby(variables, ascending=ascending)
        return self._from_temp_dataset(ds)

    def quantile(
        self,
        q: ArrayLike,
        dim: Dims = None,
        *,
        method: QuantileMethods = "linear",
        keep_attrs: bool | None = None,
        skipna: bool | None = None,
        interpolation: QuantileMethods | None = None,
    ) -> Self:
        """Compute the qth quantile of the data along the specified dimension.

        Returns the qth quantiles(s) of the array elements.

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

            The first three methods are discontiuous. The following discontinuous
            variations of the default "linear" (7.) option are also available:

                * "lower"
                * "higher"
                * "midpoint"
                * "nearest"

            See :py:func:`numpy.quantile` or [1]_ for details. The "method" argument
            was previously called "interpolation", renamed in accordance with numpy
            version 1.22.0.

        keep_attrs : bool or None, optional
            If True, the dataset's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or skipna=True has not been
            implemented (object, datetime64 or timedelta64).

        Returns
        -------
        quantiles : DataArray
            If `q` is a single quantile, then the result
            is a scalar. If multiple percentiles are given, first axis of
            the result corresponds to the quantile and a quantile dimension
            is added to the return array. The other dimensions are the
            dimensions that remain after the reduction of the array.

        See Also
        --------
        numpy.nanquantile, numpy.quantile, pandas.Series.quantile, Dataset.quantile

        Examples
        --------
        >>> da = xr.DataArray(
        ...     data=[[0.7, 4.2, 9.4, 1.5], [6.5, 7.3, 2.6, 1.9]],
        ...     coords={"x": [7, 9], "y": [1, 1.5, 2, 2.5]},
        ...     dims=("x", "y"),
        ... )
        >>> da.quantile(0)  # or da.quantile(0, dim=...)
        <xarray.DataArray ()> Size: 8B
        array(0.7)
        Coordinates:
            quantile  float64 8B 0.0
        >>> da.quantile(0, dim="x")
        <xarray.DataArray (y: 4)> Size: 32B
        array([0.7, 4.2, 2.6, 1.5])
        Coordinates:
          * y         (y) float64 32B 1.0 1.5 2.0 2.5
            quantile  float64 8B 0.0
        >>> da.quantile([0, 0.5, 1])
        <xarray.DataArray (quantile: 3)> Size: 24B
        array([0.7, 3.4, 9.4])
        Coordinates:
          * quantile  (quantile) float64 24B 0.0 0.5 1.0
        >>> da.quantile([0, 0.5, 1], dim="x")
        <xarray.DataArray (quantile: 3, y: 4)> Size: 96B
        array([[0.7 , 4.2 , 2.6 , 1.5 ],
               [3.6 , 5.75, 6.  , 1.7 ],
               [6.5 , 7.3 , 9.4 , 1.9 ]])
        Coordinates:
          * y         (y) float64 32B 1.0 1.5 2.0 2.5
          * quantile  (quantile) float64 24B 0.0 0.5 1.0

        References
        ----------
        .. [1] R. J. Hyndman and Y. Fan,
           "Sample quantiles in statistical packages,"
           The American Statistician, 50(4), pp. 361-365, 1996
        """

        ds = self._to_temp_dataset().quantile(
            q,
            dim=dim,
            keep_attrs=keep_attrs,
            method=method,
            skipna=skipna,
            interpolation=interpolation,
        )
        return self._from_temp_dataset(ds)

    def rank(
        self,
        dim: Hashable,
        *,
        pct: bool = False,
        keep_attrs: bool | None = None,
    ) -> Self:
        """Ranks the data.

        Equal values are assigned a rank that is the average of the ranks that
        would have been otherwise assigned to all of the values within that
        set.  Ranks begin at 1, not 0. If pct, computes percentage ranks.

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
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.

        Returns
        -------
        ranked : DataArray
            DataArray with the same coordinates and dtype 'float64'.

        Examples
        --------
        >>> arr = xr.DataArray([5, 6, 7], dims="x")
        >>> arr.rank("x")
        <xarray.DataArray (x: 3)> Size: 24B
        array([1., 2., 3.])
        Dimensions without coordinates: x
        """

        ds = self._to_temp_dataset().rank(dim, pct=pct, keep_attrs=keep_attrs)
        return self._from_temp_dataset(ds)

    def differentiate(
        self,
        coord: Hashable,
        edge_order: Literal[1, 2] = 1,
        datetime_unit: DatetimeUnitOptions = None,
    ) -> Self:
        """Differentiate the array with the second order accurate central
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
        datetime_unit : {"W", "D", "h", "m", "s", "ms", \
                         "us", "ns", "ps", "fs", "as", None}, optional
            Unit to compute gradient. Only valid for datetime coordinate. "Y" and "M" are not available as
            datetime_unit.

        Returns
        -------
        differentiated: DataArray

        See also
        --------
        numpy.gradient: corresponding numpy function

        Examples
        --------

        >>> da = xr.DataArray(
        ...     np.arange(12).reshape(4, 3),
        ...     dims=["x", "y"],
        ...     coords={"x": [0, 0.1, 1.1, 1.2]},
        ... )
        >>> da
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * x        (x) float64 32B 0.0 0.1 1.1 1.2
        Dimensions without coordinates: y
        >>>
        >>> da.differentiate("x")
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[30.        , 30.        , 30.        ],
               [27.54545455, 27.54545455, 27.54545455],
               [27.54545455, 27.54545455, 27.54545455],
               [30.        , 30.        , 30.        ]])
        Coordinates:
          * x        (x) float64 32B 0.0 0.1 1.1 1.2
        Dimensions without coordinates: y
        """
        ds = self._to_temp_dataset().differentiate(coord, edge_order, datetime_unit)
        return self._from_temp_dataset(ds)

    def integrate(
        self,
        coord: Hashable | Sequence[Hashable] = None,
        datetime_unit: DatetimeUnitOptions = None,
    ) -> Self:
        """Integrate along the given coordinate using the trapezoidal rule.

        .. note::
            This feature is limited to simple cartesian geometry, i.e. coord
            must be one dimensional.

        Parameters
        ----------
        coord : Hashable, or sequence of Hashable
            Coordinate(s) used for the integration.
        datetime_unit : {'W', 'D', 'h', 'm', 's', 'ms', 'us', 'ns', \
                        'ps', 'fs', 'as', None}, optional
            Specify the unit if a datetime coordinate is used.

        Returns
        -------
        integrated : DataArray

        See also
        --------
        Dataset.integrate
        numpy.trapz : corresponding numpy function

        Examples
        --------

        >>> da = xr.DataArray(
        ...     np.arange(12).reshape(4, 3),
        ...     dims=["x", "y"],
        ...     coords={"x": [0, 0.1, 1.1, 1.2]},
        ... )
        >>> da
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * x        (x) float64 32B 0.0 0.1 1.1 1.2
        Dimensions without coordinates: y
        >>>
        >>> da.integrate("x")
        <xarray.DataArray (y: 3)> Size: 24B
        array([5.4, 6.6, 7.8])
        Dimensions without coordinates: y
        """
        ds = self._to_temp_dataset().integrate(coord, datetime_unit)
        return self._from_temp_dataset(ds)

    def cumulative_integrate(
        self,
        coord: Hashable | Sequence[Hashable] = None,
        datetime_unit: DatetimeUnitOptions = None,
    ) -> Self:
        """Integrate cumulatively along the given coordinate using the trapezoidal rule.

        .. note::
            This feature is limited to simple cartesian geometry, i.e. coord
            must be one dimensional.

            The first entry of the cumulative integral is always 0, in order to keep the
            length of the dimension unchanged between input and output.

        Parameters
        ----------
        coord : Hashable, or sequence of Hashable
            Coordinate(s) used for the integration.
        datetime_unit : {'W', 'D', 'h', 'm', 's', 'ms', 'us', 'ns', \
                        'ps', 'fs', 'as', None}, optional
            Specify the unit if a datetime coordinate is used.

        Returns
        -------
        integrated : DataArray

        See also
        --------
        Dataset.cumulative_integrate
        scipy.integrate.cumulative_trapezoid : corresponding scipy function

        Examples
        --------

        >>> da = xr.DataArray(
        ...     np.arange(12).reshape(4, 3),
        ...     dims=["x", "y"],
        ...     coords={"x": [0, 0.1, 1.1, 1.2]},
        ... )
        >>> da
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[ 0,  1,  2],
               [ 3,  4,  5],
               [ 6,  7,  8],
               [ 9, 10, 11]])
        Coordinates:
          * x        (x) float64 32B 0.0 0.1 1.1 1.2
        Dimensions without coordinates: y
        >>>
        >>> da.cumulative_integrate("x")
        <xarray.DataArray (x: 4, y: 3)> Size: 96B
        array([[0.  , 0.  , 0.  ],
               [0.15, 0.25, 0.35],
               [4.65, 5.75, 6.85],
               [5.4 , 6.6 , 7.8 ]])
        Coordinates:
          * x        (x) float64 32B 0.0 0.1 1.1 1.2
        Dimensions without coordinates: y
        """
        ds = self._to_temp_dataset().cumulative_integrate(coord, datetime_unit)
        return self._from_temp_dataset(ds)

    def unify_chunks(self) -> Self:
        """Unify chunk size along all chunked dimensions of this DataArray.

        Returns
        -------
        DataArray with consistent chunk sizes for all dask-array variables

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
        Apply a function to each block of this DataArray.

        .. warning::
            This method is experimental and its signature may change.

        Parameters
        ----------
        func : callable
            User-provided function that accepts a DataArray as its first
            parameter. The function will receive a subset or 'block' of this DataArray (see below),
            corresponding to one chunk along each chunked dimension. ``func`` will be
            executed as ``func(subset_dataarray, *subset_args, **kwargs)``.

            This function must return either a single DataArray or a single Dataset.

            This function cannot add a new chunked dimension.
        args : sequence
            Passed to func after unpacking and subsetting any xarray objects by blocks.
            xarray objects in args must be aligned with this object, otherwise an error is raised.
        kwargs : mapping
            Passed verbatim to func after unpacking. xarray objects, if any, will not be
            subset to blocks. Passing dask collections in kwargs is not allowed.
        template : DataArray or Dataset, optional
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
        :func:`xarray.apply_ufunc <xarray.apply_ufunc>`
        :func:`xarray.Dataset.map_blocks <xarray.Dataset.map_blocks>`

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
        >>> array.map_blocks(calculate_anomaly, template=array).compute()
        <xarray.DataArray (time: 24)> Size: 192B
        array([ 0.12894847,  0.11323072, -0.0855964 , -0.09334032,  0.26848862,
                0.12382735,  0.22460641,  0.07650108, -0.07673453, -0.22865714,
               -0.19063865,  0.0590131 , -0.12894847, -0.11323072,  0.0855964 ,
                0.09334032, -0.26848862, -0.12382735, -0.22460641, -0.07650108,
                0.07673453,  0.22865714,  0.19063865, -0.0590131 ])
        Coordinates:
          * time     (time) object 192B 1990-01-31 00:00:00 ... 1991-12-31 00:00:00
            month    (time) int64 192B 1 2 3 4 5 6 7 8 9 10 ... 3 4 5 6 7 8 9 10 11 12

        Note that one must explicitly use ``args=[]`` and ``kwargs={}`` to pass arguments
        to the function being applied in ``xr.map_blocks()``:

        >>> array.map_blocks(
        ...     calculate_anomaly, kwargs={"groupby_type": "time.year"}, template=array
        ... )  # doctest: +ELLIPSIS
        <xarray.DataArray (time: 24)> Size: 192B
        dask.array<<this-array>-calculate_anomaly, shape=(24,), dtype=float64, chunksize=(24,), chunktype=numpy.ndarray>
        Coordinates:
          * time     (time) object 192B 1990-01-31 00:00:00 ... 1991-12-31 00:00:00
            month    (time) int64 192B dask.array<chunksize=(24,), meta=np.ndarray>
        """
        from xarray.core.parallel import map_blocks

        return map_blocks(func, self, args, kwargs, template)

    def polyfit(
        self,
        dim: Hashable,
        deg: int,
        skipna: bool | None = None,
        rcond: float | None = None,
        w: Hashable | Any | None = None,
        full: bool = False,
        cov: bool | Literal["unscaled"] = False,
    ) -> Dataset:
        """
        Least squares polynomial fit.

        This replicates the behaviour of `numpy.polyfit` but differs by skipping
        invalid values when `skipna = True`.

        Parameters
        ----------
        dim : Hashable
            Coordinate along which to fit the polynomials.
        deg : int
            Degree of the fitting polynomial.
        skipna : bool or None, optional
            If True, removes all invalid values before fitting each 1D slices of the array.
            Default is True if data is stored in a dask.array or if there is any
            invalid values, False otherwise.
        rcond : float or None, optional
            Relative condition number to the fit.
        w : Hashable, array-like or None, optional
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
            A single dataset which contains:

            polyfit_coefficients
                The coefficients of the best fit.
            polyfit_residuals
                The residuals of the least-square computation (only included if `full=True`).
                When the matrix rank is deficient, np.nan is returned.
            [dim]_matrix_rank
                The effective rank of the scaled Vandermonde coefficient matrix (only included if `full=True`)
            [dim]_singular_value
                The singular values of the scaled Vandermonde coefficient matrix (only included if `full=True`)
            polyfit_covariance
                The covariance matrix of the polynomial coefficient estimates (only included if `full=False` and `cov=True`)

        See Also
        --------
        numpy.polyfit
        numpy.polyval
        xarray.polyval
        DataArray.curvefit
        """
        # For DataArray, use the original implementation by converting to a dataset
        return self._to_temp_dataset().polyfit(
            dim, deg, skipna=skipna, rcond=rcond, w=w, full=full, cov=cov
        )

    def pad(
        self,
        pad_width: Mapping[Any, int | tuple[int, int]] | None = None,
        mode: PadModeOptions = "constant",
        stat_length: (
            int | tuple[int, int] | Mapping[Any, tuple[int, int]] | None
        ) = None,
        constant_values: (
            float | tuple[float, float] | Mapping[Any, tuple[float, float]] | None
        ) = None,
        end_values: int | tuple[int, int] | Mapping[Any, tuple[int, int]] | None = None,
        reflect_type: PadReflectOptions = None,
        keep_attrs: bool | None = None,
        **pad_width_kwargs: Any,
    ) -> Self:
        """Pad this array along one or more dimensions.

        .. warning::
            This function is experimental and its behaviour is likely to change
            especially regarding padding of dimension coordinates (or IndexVariables).

        When using one of the modes ("edge", "reflect", "symmetric", "wrap"),
        coordinates will be padded with the same mode, otherwise coordinates
        are padded using the "constant" mode with fill_value dtypes.NA.

        Parameters
        ----------
        pad_width : mapping of Hashable to tuple of int
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

        stat_length : int, tuple or mapping of Hashable to tuple, default: None
            Used in 'maximum', 'mean', 'median', and 'minimum'.  Number of
            values at edge of each axis used to calculate the statistic value.
            {dim_1: (before_1, after_1), ... dim_N: (before_N, after_N)} unique
            statistic lengths along each dimension.
            ((before, after),) yields same before and after statistic lengths
            for each dimension.
            (stat_length,) or int is a shortcut for before = after = statistic
            length for all axes.
            Default is ``None``, to use the entire axis.
        constant_values : scalar, tuple or mapping of Hashable to tuple, default: 0
            Used in 'constant'.  The values to set the padded values for each
            axis.
            ``{dim_1: (before_1, after_1), ... dim_N: (before_N, after_N)}`` unique
            pad constants along each dimension.
            ``((before, after),)`` yields same before and after constants for each
            dimension.
            ``(constant,)`` or ``constant`` is a shortcut for ``before = after = constant`` for
            all dimensions.
            Default is 0.
        end_values : scalar, tuple or mapping of Hashable to tuple, default: 0
            Used in 'linear_ramp'.  The values used for the ending value of the
            linear_ramp and that will form the edge of the padded array.
            ``{dim_1: (before_1, after_1), ... dim_N: (before_N, after_N)}`` unique
            end values along each dimension.
            ``((before, after),)`` yields same before and after end values for each
            axis.
            ``(constant,)`` or ``constant`` is a shortcut for ``before = after = constant`` for
            all axes.
            Default is 0.
        reflect_type : {"even", "odd", None}, optional
            Used in "reflect", and "symmetric". The "even" style is the
            default with an unaltered reflection around the edge value. For
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
        padded : DataArray
            DataArray with the padded coordinates and data.

        See Also
        --------
        DataArray.shift, DataArray.roll, DataArray.bfill, DataArray.ffill, numpy.pad, dask.array.pad

        Notes
        -----
        For ``mode="constant"`` and ``constant_values=None``, integer types will be
        promoted to ``float`` and padded with ``np.nan``.

        Padding coordinates will drop their corresponding index (if any) and will reset default
        indexes for dimension coordinates.

        Examples
        --------
        >>> arr = xr.DataArray([5, 6, 7], coords=[("x", [0, 1, 2])])
        >>> arr.pad(x=(1, 2), constant_values=0)
        <xarray.DataArray (x: 6)> Size: 48B
        array([0, 5, 6, 7, 0, 0])
        Coordinates:
          * x        (x) float64 48B nan 0.0 1.0 2.0 nan nan

        >>> da = xr.DataArray(
        ...     [[0, 1, 2, 3], [10, 11, 12, 13]],
        ...     dims=["x", "y"],
        ...     coords={"x": [0, 1], "y": [10, 20, 30, 40], "z": ("x", [100, 200])},
        ... )
        >>> da.pad(x=1)
        <xarray.DataArray (x: 4, y: 4)> Size: 128B
        array([[nan, nan, nan, nan],
               [ 0.,  1.,  2.,  3.],
               [10., 11., 12., 13.],
               [nan, nan, nan, nan]])
        Coordinates:
          * x        (x) float64 32B nan 0.0 1.0 nan
          * y        (y) int64 32B 10 20 30 40
            z        (x) float64 32B nan 100.0 200.0 nan

        Careful, ``constant_values`` are coerced to the data type of the array which may
        lead to a loss of precision:

        >>> da.pad(x=1, constant_values=1.23456789)
        <xarray.DataArray (x: 4, y: 4)> Size: 128B
        array([[ 1,  1,  1,  1],
               [ 0,  1,  2,  3],
               [10, 11, 12, 13],
               [ 1,  1,  1,  1]])
        Coordinates:
          * x        (x) float64 32B nan 0.0 1.0 nan
          * y        (y) int64 32B 10 20 30 40
            z        (x) float64 32B nan 100.0 200.0 nan
        """
        ds = self._to_temp_dataset().pad(
            pad_width=pad_width,
            mode=mode,
            stat_length=stat_length,
            constant_values=constant_values,
            end_values=end_values,
            reflect_type=reflect_type,
            keep_attrs=keep_attrs,
            **pad_width_kwargs,
        )
        return self._from_temp_dataset(ds)

    def idxmin(
        self,
        dim: Hashable | None = None,
        *,
        skipna: bool | None = None,
        fill_value: Any = dtypes.NA,
        keep_attrs: bool | None = None,
    ) -> Self:
        """Return the coordinate label of the minimum value along a dimension.

        Returns a new `DataArray` named after the dimension with the values of
        the coordinate labels along that dimension corresponding to minimum
        values along that dimension.

        In comparison to :py:meth:`~DataArray.argmin`, this returns the
        coordinate label while :py:meth:`~DataArray.argmin` returns the index.

        Parameters
        ----------
        dim : str, optional
            Dimension over which to apply `idxmin`.  This is optional for 1D
            arrays, but required for arrays with 2 or more dimensions.
        skipna : bool or None, default: None
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
        reduced : DataArray
            New `DataArray` object with `idxmin` applied to its data and the
            indicated dimension removed.

        See Also
        --------
        Dataset.idxmin, DataArray.idxmax, DataArray.min, DataArray.argmin

        Examples
        --------
        >>> array = xr.DataArray(
        ...     [0, 2, 1, 0, -2], dims="x", coords={"x": ["a", "b", "c", "d", "e"]}
        ... )
        >>> array.min()
        <xarray.DataArray ()> Size: 8B
        array(-2)
        >>> array.argmin(...)
        {'x': <xarray.DataArray ()> Size: 8B
        array(4)}
        >>> array.idxmin()
        <xarray.DataArray 'x' ()> Size: 4B
        array('e', dtype='<U1')

        >>> array = xr.DataArray(
        ...     [
        ...         [2.0, 1.0, 2.0, 0.0, -2.0],
        ...         [-4.0, np.nan, 2.0, np.nan, -2.0],
        ...         [np.nan, np.nan, 1.0, np.nan, np.nan],
        ...     ],
        ...     dims=["y", "x"],
        ...     coords={"y": [-1, 0, 1], "x": np.arange(5.0) ** 2},
        ... )
        >>> array.min(dim="x")
        <xarray.DataArray (y: 3)> Size: 24B
        array([-2., -4.,  1.])
        Coordinates:
          * y        (y) int64 24B -1 0 1
        >>> array.argmin(dim="x")
        <xarray.DataArray (y: 3)> Size: 24B
        array([4, 0, 2])
        Coordinates:
          * y        (y) int64 24B -1 0 1
        >>> array.idxmin(dim="x")
        <xarray.DataArray 'x' (y: 3)> Size: 24B
        array([16.,  0.,  4.])
        Coordinates:
          * y        (y) int64 24B -1 0 1
        """
        return computation._calc_idxminmax(
            array=self,
            func=lambda x, *args, **kwargs: x.argmin(*args, **kwargs),
            dim=dim,
            skipna=skipna,
            fill_value=fill_value,
            keep_attrs=keep_attrs,
        )

    def idxmax(
        self,
        dim: Hashable = None,
        *,
        skipna: bool | None = None,
        fill_value: Any = dtypes.NA,
        keep_attrs: bool | None = None,
    ) -> Self:
        """Return the coordinate label of the maximum value along a dimension.

        Returns a new `DataArray` named after the dimension with the values of
        the coordinate labels along that dimension corresponding to maximum
        values along that dimension.

        In comparison to :py:meth:`~DataArray.argmax`, this returns the
        coordinate label while :py:meth:`~DataArray.argmax` returns the index.

        Parameters
        ----------
        dim : Hashable, optional
            Dimension over which to apply `idxmax`.  This is optional for 1D
            arrays, but required for arrays with 2 or more dimensions.
        skipna : bool or None, default: None
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
        reduced : DataArray
            New `DataArray` object with `idxmax` applied to its data and the
            indicated dimension removed.

        See Also
        --------
        Dataset.idxmax, DataArray.idxmin, DataArray.max, DataArray.argmax

        Examples
        --------
        >>> array = xr.DataArray(
        ...     [0, 2, 1, 0, -2], dims="x", coords={"x": ["a", "b", "c", "d", "e"]}
        ... )
        >>> array.max()
        <xarray.DataArray ()> Size: 8B
        array(2)
        >>> array.argmax(...)
        {'x': <xarray.DataArray ()> Size: 8B
        array(1)}
        >>> array.idxmax()
        <xarray.DataArray 'x' ()> Size: 4B
        array('b', dtype='<U1')

        >>> array = xr.DataArray(
        ...     [
        ...         [2.0, 1.0, 2.0, 0.0, -2.0],
        ...         [-4.0, np.nan, 2.0, np.nan, -2.0],
        ...         [np.nan, np.nan, 1.0, np.nan, np.nan],
        ...     ],
        ...     dims=["y", "x"],
        ...     coords={"y": [-1, 0, 1], "x": np.arange(5.0) ** 2},
        ... )
        >>> array.max(dim="x")
        <xarray.DataArray (y: 3)> Size: 24B
        array([2., 2., 1.])
        Coordinates:
          * y        (y) int64 24B -1 0 1
        >>> array.argmax(dim="x")
        <xarray.DataArray (y: 3)> Size: 24B
        array([0, 2, 2])
        Coordinates:
          * y        (y) int64 24B -1 0 1
        >>> array.idxmax(dim="x")
        <xarray.DataArray 'x' (y: 3)> Size: 24B
        array([0., 4., 4.])
        Coordinates:
          * y        (y) int64 24B -1 0 1
        """
        return computation._calc_idxminmax(
            array=self,
            func=lambda x, *args, **kwargs: x.argmax(*args, **kwargs),
            dim=dim,
            skipna=skipna,
            fill_value=fill_value,
            keep_attrs=keep_attrs,
        )

    def argmin(
        self,
        dim: Dims = None,
        *,
        axis: int | None = None,
        keep_attrs: bool | None = None,
        skipna: bool | None = None,
    ) -> Self | dict[Hashable, Self]:
        """Index or indices of the minimum of the DataArray over one or more dimensions.

        If a sequence is passed to 'dim', then result returned as dict of DataArrays,
        which can be passed directly to isel(). If a single str is passed to 'dim' then
        returns a DataArray with dtype int.

        If there are multiple minima, the indices of the first one found will be
        returned.

        Parameters
        ----------
        dim : "...", str, Iterable of Hashable or None, optional
            The dimensions over which to find the minimum. By default, finds minimum over
            all dimensions - for now returning an int for backward compatibility, but
            this is deprecated, in future will return a dict with indices for all
            dimensions; to return a dict with all dimensions now, pass '...'.
        axis : int or None, optional
            Axis over which to apply `argmin`. Only one of the 'dim' and 'axis' arguments
            can be supplied.
        keep_attrs : bool or None, optional
            If True, the attributes (`attrs`) will be copied from the original
            object to the new one. If False, the new object will be
            returned without attributes.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or skipna=True has not been
            implemented (object, datetime64 or timedelta64).

        Returns
        -------
        result : DataArray or dict of DataArray

        See Also
        --------
        Variable.argmin, DataArray.idxmin

        Examples
        --------
        >>> array = xr.DataArray([0, 2, -1, 3], dims="x")
        >>> array.min()
        <xarray.DataArray ()> Size: 8B
        array(-1)
        >>> array.argmin(...)
        {'x': <xarray.DataArray ()> Size: 8B
        array(2)}
        >>> array.isel(array.argmin(...))
        <xarray.DataArray ()> Size: 8B
        array(-1)

        >>> array = xr.DataArray(
        ...     [[[3, 2, 1], [3, 1, 2], [2, 1, 3]], [[1, 3, 2], [2, -5, 1], [2, 3, 1]]],
        ...     dims=("x", "y", "z"),
        ... )
        >>> array.min(dim="x")
        <xarray.DataArray (y: 3, z: 3)> Size: 72B
        array([[ 1,  2,  1],
               [ 2, -5,  1],
               [ 2,  1,  1]])
        Dimensions without coordinates: y, z
        >>> array.argmin(dim="x")
        <xarray.DataArray (y: 3, z: 3)> Size: 72B
        array([[1, 0, 0],
               [1, 1, 1],
               [0, 0, 1]])
        Dimensions without coordinates: y, z
        >>> array.argmin(dim=["x"])
        {'x': <xarray.DataArray (y: 3, z: 3)> Size: 72B
        array([[1, 0, 0],
               [1, 1, 1],
               [0, 0, 1]])
        Dimensions without coordinates: y, z}
        >>> array.min(dim=("x", "z"))
        <xarray.DataArray (y: 3)> Size: 24B
        array([ 1, -5,  1])
        Dimensions without coordinates: y
        >>> array.argmin(dim=["x", "z"])
        {'x': <xarray.DataArray (y: 3)> Size: 24B
        array([0, 1, 0])
        Dimensions without coordinates: y, 'z': <xarray.DataArray (y: 3)> Size: 24B
        array([2, 1, 1])
        Dimensions without coordinates: y}
        >>> array.isel(array.argmin(dim=["x", "z"]))
        <xarray.DataArray (y: 3)> Size: 24B
        array([ 1, -5,  1])
        Dimensions without coordinates: y
        """
        result = self.variable.argmin(dim, axis, keep_attrs, skipna)
        if isinstance(result, dict):
            return {k: self._replace_maybe_drop_dims(v) for k, v in result.items()}
        else:
            return self._replace_maybe_drop_dims(result)

    def argmax(
        self,
        dim: Dims = None,
        *,
        axis: int | None = None,
        keep_attrs: bool | None = None,
        skipna: bool | None = None,
    ) -> Self | dict[Hashable, Self]:
        """Index or indices of the maximum of the DataArray over one or more dimensions.

        If a sequence is passed to 'dim', then result returned as dict of DataArrays,
        which can be passed directly to isel(). If a single str is passed to 'dim' then
        returns a DataArray with dtype int.

        If there are multiple maxima, the indices of the first one found will be
        returned.

        Parameters
        ----------
        dim : "...", str, Iterable of Hashable or None, optional
            The dimensions over which to find the maximum. By default, finds maximum over
            all dimensions - for now returning an int for backward compatibility, but
            this is deprecated, in future will return a dict with indices for all
            dimensions; to return a dict with all dimensions now, pass '...'.
        axis : int or None, optional
            Axis over which to apply `argmax`. Only one of the 'dim' and 'axis' arguments
            can be supplied.
        keep_attrs : bool or None, optional
            If True, the attributes (`attrs`) will be copied from the original
            object to the new one. If False, the new object will be
            returned without attributes.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or skipna=True has not been
            implemented (object, datetime64 or timedelta64).

        Returns
        -------
        result : DataArray or dict of DataArray

        See Also
        --------
        Variable.argmax, DataArray.idxmax

        Examples
        --------
        >>> array = xr.DataArray([0, 2, -1, 3], dims="x")
        >>> array.max()
        <xarray.DataArray ()> Size: 8B
        array(3)
        >>> array.argmax(...)
        {'x': <xarray.DataArray ()> Size: 8B
        array(3)}
        >>> array.isel(array.argmax(...))
        <xarray.DataArray ()> Size: 8B
        array(3)

        >>> array = xr.DataArray(
        ...     [[[3, 2, 1], [3, 1, 2], [2, 1, 3]], [[1, 3, 2], [2, 5, 1], [2, 3, 1]]],
        ...     dims=("x", "y", "z"),
        ... )
        >>> array.max(dim="x")
        <xarray.DataArray (y: 3, z: 3)> Size: 72B
        array([[3, 3, 2],
               [3, 5, 2],
               [2, 3, 3]])
        Dimensions without coordinates: y, z
        >>> array.argmax(dim="x")
        <xarray.DataArray (y: 3, z: 3)> Size: 72B
        array([[0, 1, 1],
               [0, 1, 0],
               [0, 1, 0]])
        Dimensions without coordinates: y, z
        >>> array.argmax(dim=["x"])
        {'x': <xarray.DataArray (y: 3, z: 3)> Size: 72B
        array([[0, 1, 1],
               [0, 1, 0],
               [0, 1, 0]])
        Dimensions without coordinates: y, z}
        >>> array.max(dim=("x", "z"))
        <xarray.DataArray (y: 3)> Size: 24B
        array([3, 5, 3])
        Dimensions without coordinates: y
        >>> array.argmax(dim=["x", "z"])
        {'x': <xarray.DataArray (y: 3)> Size: 24B
        array([0, 1, 0])
        Dimensions without coordinates: y, 'z': <xarray.DataArray (y: 3)> Size: 24B
        array([0, 1, 2])
        Dimensions without coordinates: y}
        >>> array.isel(array.argmax(dim=["x", "z"]))
        <xarray.DataArray (y: 3)> Size: 24B
        array([3, 5, 3])
        Dimensions without coordinates: y
        """
        result = self.variable.argmax(dim, axis, keep_attrs, skipna)
        if isinstance(result, dict):
            return {k: self._replace_maybe_drop_dims(v) for k, v in result.items()}
        else:
            return self._replace_maybe_drop_dims(result)

    def query(
        self,
        queries: Mapping[Any, Any] | None = None,
        parser: QueryParserOptions = "pandas",
        engine: QueryEngineOptions = None,
        missing_dims: ErrorOptionsWithWarn = "raise",
        **queries_kwargs: Any,
    ) -> DataArray:
        """Return a new data array indexed along the specified
        dimension(s), where the indexers are given as strings containing
        Python expressions to be evaluated against the values in the array.

        Parameters
        ----------
        queries : dict-like or None, optional
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
            DataArray:

            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        **queries_kwargs : {dim: query, ...}, optional
            The keyword arguments form of ``queries``.
            One of queries or queries_kwargs must be provided.

        Returns
        -------
        obj : DataArray
            A new DataArray with the same contents as this dataset, indexed by
            the results of the appropriate queries.

        See Also
        --------
        DataArray.isel
        Dataset.query
        pandas.eval

        Examples
        --------
        >>> da = xr.DataArray(np.arange(0, 5, 1), dims="x", name="a")
        >>> da
        <xarray.DataArray 'a' (x: 5)> Size: 40B
        array([0, 1, 2, 3, 4])
        Dimensions without coordinates: x
        >>> da.query(x="a > 2")
        <xarray.DataArray 'a' (x: 2)> Size: 16B
        array([3, 4])
        Dimensions without coordinates: x
        """

        ds = self._to_dataset_whole(shallow_copy=True)
        ds = ds.query(
            queries=queries,
            parser=parser,
            engine=engine,
            missing_dims=missing_dims,
            **queries_kwargs,
        )
        return ds[self.name]

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
    ) -> Dataset:
        """
        Curve fitting optimization for arbitrary functions.

        Wraps :py:func:`scipy.optimize.curve_fit` with :py:func:`~xarray.apply_ufunc`.

        Parameters
        ----------
        coords : Hashable, DataArray, or sequence of DataArray or Hashable
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
        p0 : dict-like or None, optional
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
        param_names : sequence of Hashable or None, optional
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

        Examples
        --------
        Generate some exponentially decaying data, where the decay constant and amplitude are
        different for different values of the coordinate ``x``:

        >>> rng = np.random.default_rng(seed=0)
        >>> def exp_decay(t, time_constant, amplitude):
        ...     return np.exp(-t / time_constant) * amplitude
        ...
        >>> t = np.arange(11)
        >>> da = xr.DataArray(
        ...     np.stack(
        ...         [
        ...             exp_decay(t, 1, 0.1),
        ...             exp_decay(t, 2, 0.2),
        ...             exp_decay(t, 3, 0.3),
        ...         ]
        ...     )
        ...     + rng.normal(size=(3, t.size)) * 0.01,
        ...     coords={"x": [0, 1, 2], "time": t},
        ... )
        >>> da
        <xarray.DataArray (x: 3, time: 11)> Size: 264B
        array([[ 0.1012573 ,  0.0354669 ,  0.01993775,  0.00602771, -0.00352513,
                 0.00428975,  0.01328788,  0.009562  , -0.00700381, -0.01264187,
                -0.0062282 ],
               [ 0.20041326,  0.09805582,  0.07138797,  0.03216692,  0.01974438,
                 0.01097441,  0.00679441,  0.01015578,  0.01408826,  0.00093645,
                 0.01501222],
               [ 0.29334805,  0.21847449,  0.16305984,  0.11130396,  0.07164415,
                 0.04744543,  0.03602333,  0.03129354,  0.01074885,  0.01284436,
                 0.00910995]])
        Coordinates:
          * x        (x) int64 24B 0 1 2
          * time     (time) int64 88B 0 1 2 3 4 5 6 7 8 9 10

        Fit the exponential decay function to the data along the ``time`` dimension:

        >>> fit_result = da.curvefit("time", exp_decay)
        >>> fit_result["curvefit_coefficients"].sel(
        ...     param="time_constant"
        ... )  # doctest: +NUMBER
        <xarray.DataArray 'curvefit_coefficients' (x: 3)> Size: 24B
        array([1.05692036, 1.73549638, 2.94215771])
        Coordinates:
          * x        (x) int64 24B 0 1 2
            param    <U13 52B 'time_constant'
        >>> fit_result["curvefit_coefficients"].sel(param="amplitude")
        <xarray.DataArray 'curvefit_coefficients' (x: 3)> Size: 24B
        array([0.1005489 , 0.19631423, 0.30003579])
        Coordinates:
          * x        (x) int64 24B 0 1 2
            param    <U13 52B 'amplitude'

        An initial guess can also be given with the ``p0`` arg (although it does not make much
        of a difference in this simple example). To have a different guess for different
        coordinate points, the guess can be a DataArray. Here we use the same initial guess
        for the amplitude but different guesses for the time constant:

        >>> fit_result = da.curvefit(
        ...     "time",
        ...     exp_decay,
        ...     p0={
        ...         "amplitude": 0.2,
        ...         "time_constant": xr.DataArray([1, 2, 3], coords=[da.x]),
        ...     },
        ... )
        >>> fit_result["curvefit_coefficients"].sel(param="time_constant")
        <xarray.DataArray 'curvefit_coefficients' (x: 3)> Size: 24B
        array([1.0569213 , 1.73550052, 2.94215733])
        Coordinates:
          * x        (x) int64 24B 0 1 2
            param    <U13 52B 'time_constant'
        >>> fit_result["curvefit_coefficients"].sel(param="amplitude")
        <xarray.DataArray 'curvefit_coefficients' (x: 3)> Size: 24B
        array([0.10054889, 0.1963141 , 0.3000358 ])
        Coordinates:
          * x        (x) int64 24B 0 1 2
            param    <U13 52B 'amplitude'

        See Also
        --------
        DataArray.polyfit
        scipy.optimize.curve_fit
        xarray.DataArray.xlm.modelfit
            External method from `xarray-lmfit <https://xarray-lmfit.readthedocs.io/>`_
            with more curve fitting functionality.
        """
        # For DataArray, use the original implementation by converting to a dataset first
        return self._to_temp_dataset().curvefit(
            coords,
            func,
            reduce_dims=reduce_dims,
            skipna=skipna,
            p0=p0,
            bounds=bounds,
            param_names=param_names,
            errors=errors,
            kwargs=kwargs,
        )

    def drop_duplicates(
        self,
        dim: Hashable | Iterable[Hashable],
        *,
        keep: Literal["first", "last", False] = "first",
    ) -> Self:
        """Returns a new DataArray with duplicate dimension values removed.

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
        DataArray

        See Also
        --------
        Dataset.drop_duplicates

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.arange(25).reshape(5, 5),
        ...     dims=("x", "y"),
        ...     coords={"x": np.array([0, 0, 1, 2, 3]), "y": np.array([0, 1, 2, 3, 3])},
        ... )
        >>> da
        <xarray.DataArray (x: 5, y: 5)> Size: 200B
        array([[ 0,  1,  2,  3,  4],
               [ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Coordinates:
          * x        (x) int64 40B 0 0 1 2 3
          * y        (y) int64 40B 0 1 2 3 3

        >>> da.drop_duplicates(dim="x")
        <xarray.DataArray (x: 4, y: 5)> Size: 160B
        array([[ 0,  1,  2,  3,  4],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
          * y        (y) int64 40B 0 1 2 3 3

        >>> da.drop_duplicates(dim="x", keep="last")
        <xarray.DataArray (x: 4, y: 5)> Size: 160B
        array([[ 5,  6,  7,  8,  9],
               [10, 11, 12, 13, 14],
               [15, 16, 17, 18, 19],
               [20, 21, 22, 23, 24]])
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
          * y        (y) int64 40B 0 1 2 3 3

        Drop all duplicate dimension values:

        >>> da.drop_duplicates(dim=...)
        <xarray.DataArray (x: 4, y: 4)> Size: 128B
        array([[ 0,  1,  2,  3],
               [10, 11, 12, 13],
               [15, 16, 17, 18],
               [20, 21, 22, 23]])
        Coordinates:
          * x        (x) int64 32B 0 1 2 3
          * y        (y) int64 32B 0 1 2 3
        """
        deduplicated = self._to_temp_dataset().drop_duplicates(dim, keep=keep)
        return self._from_temp_dataset(deduplicated)

    def convert_calendar(
        self,
        calendar: str,
        dim: str = "time",
        align_on: str | None = None,
        missing: Any | None = None,
        use_cftime: bool | None = None,
    ) -> Self:
        """Convert the DataArray to another calendar.

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
        dim : str
            Name of the time coordinate.
        align_on : {None, 'date', 'year'}
            Must be specified when either source or target is a `360_day` calendar,
           ignored otherwise. See Notes.
        missing : Optional[any]
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
        use_cftime : boolean, optional
            Whether to use cftime objects in the output, only used if `calendar`
            is one of {"proleptic_gregorian", "gregorian" or "standard"}.
            If True, the new time axis uses cftime objects.
            If None (default), it uses :py:class:`numpy.datetime64` values if the
            date range permits it, and :py:class:`cftime.datetime` objects if not.
            If False, it uses :py:class:`numpy.datetime64`  or fails.

        Returns
        -------
        DataArray
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
        dim: str = "time",
    ) -> Self:
        """Interpolates the DataArray to another calendar based on decimal year measure.

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
        dim : str
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
    ) -> DataArrayGroupBy:
        """Returns a DataArrayGroupBy object for performing grouped operations.

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
        eagerly_compute_group: bool, optional
            This argument is deprecated.
        **groupers : Mapping of str to Grouper or Resampler
            Mapping of variable name to group by to :py:class:`Grouper` or :py:class:`Resampler` object.
            One of ``group`` or ``groupers`` must be provided.
            Only a single ``grouper`` is allowed at present.

        Returns
        -------
        grouped : DataArrayGroupBy
            A `DataArrayGroupBy` object patterned after `pandas.GroupBy` that can be
            iterated over in the form of `(unique_value, grouped_array)` pairs.

        Examples
        --------
        Calculate daily anomalies for daily data:

        >>> da = xr.DataArray(
        ...     np.linspace(0, 1826, num=1827),
        ...     coords=[pd.date_range("2000-01-01", "2004-12-31", freq="D")],
        ...     dims="time",
        ... )
        >>> da
        <xarray.DataArray (time: 1827)> Size: 15kB
        array([0.000e+00, 1.000e+00, 2.000e+00, ..., 1.824e+03, 1.825e+03,
               1.826e+03], shape=(1827,))
        Coordinates:
          * time     (time) datetime64[ns] 15kB 2000-01-01 2000-01-02 ... 2004-12-31
        >>> da.groupby("time.dayofyear") - da.groupby("time.dayofyear").mean("time")
        <xarray.DataArray (time: 1827)> Size: 15kB
        array([-730.8, -730.8, -730.8, ...,  730.2,  730.2,  730.5], shape=(1827,))
        Coordinates:
          * time       (time) datetime64[ns] 15kB 2000-01-01 2000-01-02 ... 2004-12-31
            dayofyear  (time) int64 15kB 1 2 3 4 5 6 7 8 ... 360 361 362 363 364 365 366

        Use a ``Grouper`` object to be more explicit

        >>> da.coords["dayofyear"] = da.time.dt.dayofyear
        >>> da.groupby(dayofyear=xr.groupers.UniqueGrouper()).mean()
        <xarray.DataArray (dayofyear: 366)> Size: 3kB
        array([ 730.8,  731.8,  732.8, ..., 1093.8, 1094.8, 1095.5])
        Coordinates:
          * dayofyear  (dayofyear) int64 3kB 1 2 3 4 5 6 7 ... 361 362 363 364 365 366

        >>> da = xr.DataArray(
        ...     data=np.arange(12).reshape((4, 3)),
        ...     dims=("x", "y"),
        ...     coords={"x": [10, 20, 30, 40], "letters": ("x", list("abba"))},
        ... )

        Grouping by a single variable is easy

        >>> da.groupby("letters")
        <DataArrayGroupBy, grouped over 1 grouper(s), 2 groups in total:
            'letters': UniqueGrouper('letters'), 2/2 groups with labels 'a', 'b'>

        Execute a reduction

        >>> da.groupby("letters").sum()
        <xarray.DataArray (letters: 2, y: 3)> Size: 48B
        array([[ 9, 11, 13],
               [ 9, 11, 13]])
        Coordinates:
          * letters  (letters) object 16B 'a' 'b'
        Dimensions without coordinates: y

        Grouping by multiple variables

        >>> da.groupby(["letters", "x"])
        <DataArrayGroupBy, grouped over 2 grouper(s), 8 groups in total:
            'letters': UniqueGrouper('letters'), 2/2 groups with labels 'a', 'b'
            'x': UniqueGrouper('x'), 4/4 groups with labels 10, 20, 30, 40>

        Use Grouper objects to express more complicated GroupBy operations

        >>> from xarray.groupers import BinGrouper, UniqueGrouper
        >>>
        >>> da.groupby(x=BinGrouper(bins=[5, 15, 25]), letters=UniqueGrouper()).sum()
        <xarray.DataArray (x_bins: 2, letters: 2, y: 3)> Size: 96B
        array([[[ 0.,  1.,  2.],
                [nan, nan, nan]],
        <BLANKLINE>
               [[nan, nan, nan],
                [ 3.,  4.,  5.]]])
        Coordinates:
          * x_bins   (x_bins) interval[int64, right] 32B (5, 15] (15, 25]
          * letters  (letters) object 16B 'a' 'b'
        Dimensions without coordinates: y


        See Also
        --------
        :ref:`groupby`
            Users guide explanation of how to group and bin data.

        :doc:`xarray-tutorial:intermediate/computation/01-high-level-computation-patterns`
            Tutorial on :py:func:`~xarray.DataArray.Groupby` for windowed computation

        :doc:`xarray-tutorial:fundamentals/03.2_groupby_with_xarray`
            Tutorial on :py:func:`~xarray.DataArray.Groupby` demonstrating reductions, transformation and comparison with :py:func:`~xarray.DataArray.resample`

        :external:py:meth:`pandas.DataFrame.groupby <pandas.DataFrame.groupby>`
        :func:`DataArray.groupby_bins <DataArray.groupby_bins>`
        :func:`Dataset.groupby <Dataset.groupby>`
        :func:`core.groupby.DataArrayGroupBy <core.groupby.DataArrayGroupBy>`
        :func:`DataArray.coarsen <DataArray.coarsen>`
        :func:`Dataset.resample <Dataset.resample>`
        :func:`DataArray.resample <DataArray.resample>`
        """
        from xarray.core.groupby import (
            DataArrayGroupBy,
            _parse_group_and_groupers,
            _validate_groupby_squeeze,
        )

        _validate_groupby_squeeze(squeeze)
        rgroupers = _parse_group_and_groupers(
            self, group, groupers, eagerly_compute_group=eagerly_compute_group
        )
        return DataArrayGroupBy(self, rgroupers, restore_coord_dims=restore_coord_dims)

    @_deprecate_positional_args("v2024.07.0")
    def groupby_bins(
        self,
        group: Hashable | DataArray | IndexVariable,
        bins: Bins,
        right: bool = True,
        labels: ArrayLike | Literal[False] | None = None,
        precision: int = 3,
        include_lowest: bool = False,
        squeeze: Literal[False] = False,
        restore_coord_dims: bool = False,
        duplicates: Literal["raise", "drop"] = "raise",
        eagerly_compute_group: Literal[False] | None = None,
    ) -> DataArrayGroupBy:
        """Returns a DataArrayGroupBy object for performing grouped operations.

        Rather than using all unique values of `group`, the values are discretized
        first by applying `pandas.cut` [1]_ to `group`.

        Parameters
        ----------
        group : Hashable, DataArray or IndexVariable
            Array whose binned values should be used to group this array. If a
            Hashable, must be the name of a coordinate contained in this dataarray.
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
        labels : array-like, False or None, default: None
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
        eagerly_compute_group: bool, optional
            This argument is deprecated.

        Returns
        -------
        grouped : DataArrayGroupBy
            A `DataArrayGroupBy` object patterned after `pandas.GroupBy` that can be
            iterated over in the form of `(unique_value, grouped_array)` pairs.
            The name of the group has the added suffix `_bins` in order to
            distinguish it from the original variable.

        See Also
        --------
        :ref:`groupby`
            Users guide explanation of how to group and bin data.
        DataArray.groupby
        Dataset.groupby_bins
        core.groupby.DataArrayGroupBy
        pandas.DataFrame.groupby

        References
        ----------
        .. [1] https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.cut.html
        """
        from xarray.core.groupby import (
            DataArrayGroupBy,
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

        return DataArrayGroupBy(
            self,
            (rgrouper,),
            restore_coord_dims=restore_coord_dims,
        )

    def weighted(self, weights: DataArray) -> DataArrayWeighted:
        """
        Weighted DataArray operations.

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
        computation.weighted.DataArrayWeighted

        See Also
        --------
        :func:`Dataset.weighted <Dataset.weighted>`

        :ref:`compute.weighted`
            User guide on weighted array reduction using :py:func:`~xarray.DataArray.weighted`

        :doc:`xarray-tutorial:fundamentals/03.4_weighted`
            Tutorial on Weighted Reduction using :py:func:`~xarray.DataArray.weighted`

        """
        from xarray.computation.weighted import DataArrayWeighted

        return DataArrayWeighted(self, weights)

    def rolling(
        self,
        dim: Mapping[Any, int] | None = None,
        min_periods: int | None = None,
        center: bool | Mapping[Any, bool] = False,
        **window_kwargs: int,
    ) -> DataArrayRolling:
        """
        Rolling window object for DataArrays.

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
        computation.rolling.DataArrayRolling

        Examples
        --------
        Create rolling seasonal average of monthly data e.g. DJF, JFM, ..., SON:

        >>> da = xr.DataArray(
        ...     np.linspace(0, 11, num=12),
        ...     coords=[
        ...         pd.date_range(
        ...             "1999-12-15",
        ...             periods=12,
        ...             freq=pd.DateOffset(months=1),
        ...         )
        ...     ],
        ...     dims="time",
        ... )
        >>> da
        <xarray.DataArray (time: 12)> Size: 96B
        array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11.])
        Coordinates:
          * time     (time) datetime64[ns] 96B 1999-12-15 2000-01-15 ... 2000-11-15
        >>> da.rolling(time=3, center=True).mean()
        <xarray.DataArray (time: 12)> Size: 96B
        array([nan,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., nan])
        Coordinates:
          * time     (time) datetime64[ns] 96B 1999-12-15 2000-01-15 ... 2000-11-15

        Remove the NaNs using ``dropna()``:

        >>> da.rolling(time=3, center=True).mean().dropna("time")
        <xarray.DataArray (time: 10)> Size: 80B
        array([ 1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10.])
        Coordinates:
          * time     (time) datetime64[ns] 80B 2000-01-15 2000-02-15 ... 2000-10-15

        See Also
        --------
        DataArray.cumulative
        Dataset.rolling
        computation.rolling.DataArrayRolling
        """
        from xarray.computation.rolling import DataArrayRolling

        dim = either_dict_or_kwargs(dim, window_kwargs, "rolling")
        return DataArrayRolling(self, dim, min_periods=min_periods, center=center)

    def cumulative(
        self,
        dim: str | Iterable[Hashable],
        min_periods: int = 1,
    ) -> DataArrayRolling:
        """
        Accumulating object for DataArrays.

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
        computation.rolling.DataArrayRolling

        Examples
        --------
        Create rolling seasonal average of monthly data e.g. DJF, JFM, ..., SON:

        >>> da = xr.DataArray(
        ...     np.linspace(0, 11, num=12),
        ...     coords=[
        ...         pd.date_range(
        ...             "1999-12-15",
        ...             periods=12,
        ...             freq=pd.DateOffset(months=1),
        ...         )
        ...     ],
        ...     dims="time",
        ... )

        >>> da
        <xarray.DataArray (time: 12)> Size: 96B
        array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11.])
        Coordinates:
          * time     (time) datetime64[ns] 96B 1999-12-15 2000-01-15 ... 2000-11-15

        >>> da.cumulative("time").sum()
        <xarray.DataArray (time: 12)> Size: 96B
        array([ 0.,  1.,  3.,  6., 10., 15., 21., 28., 36., 45., 55., 66.])
        Coordinates:
          * time     (time) datetime64[ns] 96B 1999-12-15 2000-01-15 ... 2000-11-15

        See Also
        --------
        DataArray.rolling
        Dataset.cumulative
        computation.rolling.DataArrayRolling
        """
        from xarray.computation.rolling import DataArrayRolling

        # Could we abstract this "normalize and check 'dim'" logic? It's currently shared
        # with the same method in Dataset.
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

        return DataArrayRolling(self, dim, min_periods=min_periods, center=False)

    def coarsen(
        self,
        dim: Mapping[Any, int] | None = None,
        boundary: CoarsenBoundaryOptions = "exact",
        side: SideOptions | Mapping[Any, SideOptions] = "left",
        coord_func: str | Callable | Mapping[Any, str | Callable] = "mean",
        **window_kwargs: int,
    ) -> DataArrayCoarsen:
        """
        Coarsen object for DataArrays.

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
        computation.rolling.DataArrayCoarsen

        Examples
        --------
        Coarsen the long time series by averaging over every three days.

        >>> da = xr.DataArray(
        ...     np.linspace(0, 364, num=364),
        ...     dims="time",
        ...     coords={"time": pd.date_range("1999-12-15", periods=364)},
        ... )
        >>> da  # +doctest: ELLIPSIS
        <xarray.DataArray (time: 364)> Size: 3kB
        array([  0.        ,   1.00275482,   2.00550964,   3.00826446,
                 4.01101928,   5.0137741 ,   6.01652893,   7.01928375,
                 8.02203857,   9.02479339,  10.02754821,  11.03030303,
                12.03305785,  13.03581267,  14.03856749,  15.04132231,
                16.04407713,  17.04683196,  18.04958678,  19.0523416 ,
                20.05509642,  21.05785124,  22.06060606,  23.06336088,
                24.0661157 ,  25.06887052,  26.07162534,  27.07438017,
                28.07713499,  29.07988981,  30.08264463,  31.08539945,
                32.08815427,  33.09090909,  34.09366391,  35.09641873,
                36.09917355,  37.10192837,  38.1046832 ,  39.10743802,
                40.11019284,  41.11294766,  42.11570248,  43.1184573 ,
                44.12121212,  45.12396694,  46.12672176,  47.12947658,
                48.1322314 ,  49.13498623,  50.13774105,  51.14049587,
                52.14325069,  53.14600551,  54.14876033,  55.15151515,
                56.15426997,  57.15702479,  58.15977961,  59.16253444,
                60.16528926,  61.16804408,  62.1707989 ,  63.17355372,
                64.17630854,  65.17906336,  66.18181818,  67.184573  ,
                68.18732782,  69.19008264,  70.19283747,  71.19559229,
                72.19834711,  73.20110193,  74.20385675,  75.20661157,
                76.20936639,  77.21212121,  78.21487603,  79.21763085,
        ...
               284.78236915, 285.78512397, 286.78787879, 287.79063361,
               288.79338843, 289.79614325, 290.79889807, 291.80165289,
               292.80440771, 293.80716253, 294.80991736, 295.81267218,
               296.815427  , 297.81818182, 298.82093664, 299.82369146,
               300.82644628, 301.8292011 , 302.83195592, 303.83471074,
               304.83746556, 305.84022039, 306.84297521, 307.84573003,
               308.84848485, 309.85123967, 310.85399449, 311.85674931,
               312.85950413, 313.86225895, 314.86501377, 315.8677686 ,
               316.87052342, 317.87327824, 318.87603306, 319.87878788,
               320.8815427 , 321.88429752, 322.88705234, 323.88980716,
               324.89256198, 325.8953168 , 326.89807163, 327.90082645,
               328.90358127, 329.90633609, 330.90909091, 331.91184573,
               332.91460055, 333.91735537, 334.92011019, 335.92286501,
               336.92561983, 337.92837466, 338.93112948, 339.9338843 ,
               340.93663912, 341.93939394, 342.94214876, 343.94490358,
               344.9476584 , 345.95041322, 346.95316804, 347.95592287,
               348.95867769, 349.96143251, 350.96418733, 351.96694215,
               352.96969697, 353.97245179, 354.97520661, 355.97796143,
               356.98071625, 357.98347107, 358.9862259 , 359.98898072,
               360.99173554, 361.99449036, 362.99724518, 364.        ])
        Coordinates:
          * time     (time) datetime64[ns] 3kB 1999-12-15 1999-12-16 ... 2000-12-12
        >>> da.coarsen(time=3, boundary="trim").mean()  # +doctest: ELLIPSIS
        <xarray.DataArray (time: 121)> Size: 968B
        array([  1.00275482,   4.01101928,   7.01928375,  10.02754821,
                13.03581267,  16.04407713,  19.0523416 ,  22.06060606,
                25.06887052,  28.07713499,  31.08539945,  34.09366391,
                37.10192837,  40.11019284,  43.1184573 ,  46.12672176,
                49.13498623,  52.14325069,  55.15151515,  58.15977961,
                61.16804408,  64.17630854,  67.184573  ,  70.19283747,
                73.20110193,  76.20936639,  79.21763085,  82.22589532,
                85.23415978,  88.24242424,  91.25068871,  94.25895317,
                97.26721763, 100.27548209, 103.28374656, 106.29201102,
               109.30027548, 112.30853994, 115.31680441, 118.32506887,
               121.33333333, 124.3415978 , 127.34986226, 130.35812672,
               133.36639118, 136.37465565, 139.38292011, 142.39118457,
               145.39944904, 148.4077135 , 151.41597796, 154.42424242,
               157.43250689, 160.44077135, 163.44903581, 166.45730028,
               169.46556474, 172.4738292 , 175.48209366, 178.49035813,
               181.49862259, 184.50688705, 187.51515152, 190.52341598,
               193.53168044, 196.5399449 , 199.54820937, 202.55647383,
               205.56473829, 208.57300275, 211.58126722, 214.58953168,
               217.59779614, 220.60606061, 223.61432507, 226.62258953,
               229.63085399, 232.63911846, 235.64738292, 238.65564738,
               241.66391185, 244.67217631, 247.68044077, 250.68870523,
               253.6969697 , 256.70523416, 259.71349862, 262.72176309,
               265.73002755, 268.73829201, 271.74655647, 274.75482094,
               277.7630854 , 280.77134986, 283.77961433, 286.78787879,
               289.79614325, 292.80440771, 295.81267218, 298.82093664,
               301.8292011 , 304.83746556, 307.84573003, 310.85399449,
               313.86225895, 316.87052342, 319.87878788, 322.88705234,
               325.8953168 , 328.90358127, 331.91184573, 334.92011019,
               337.92837466, 340.93663912, 343.94490358, 346.95316804,
               349.96143251, 352.96969697, 355.97796143, 358.9862259 ,
               361.99449036])
        Coordinates:
          * time     (time) datetime64[ns] 968B 1999-12-16 1999-12-19 ... 2000-12-10
        >>>

        See Also
        --------
        :class:`computation.rolling.DataArrayCoarsen <computation.rolling.DataArrayCoarsen>`
        :func:`Dataset.coarsen <Dataset.coarsen>`

        :ref:`reshape.coarsen`
            User guide describing :py:func:`~xarray.DataArray.coarsen`

        :ref:`compute.coarsen`
            User guide on block aggregation :py:func:`~xarray.DataArray.coarsen`

        :doc:`xarray-tutorial:fundamentals/03.3_windowed`
            Tutorial on windowed computation using :py:func:`~xarray.DataArray.coarsen`

        """
        from xarray.computation.rolling import DataArrayCoarsen

        dim = either_dict_or_kwargs(dim, window_kwargs, "coarsen")
        return DataArrayCoarsen(
            self,
            dim,
            boundary=boundary,
            side=side,
            coord_func=coord_func,
        )

    @_deprecate_positional_args("v2024.07.0")
    def resample(
        self,
        indexer: Mapping[Hashable, ResampleCompatible | Resampler] | None = None,
        *,
        skipna: bool | None = None,
        closed: SideOptions | None = None,
        label: SideOptions | None = None,
        offset: pd.Timedelta | datetime.timedelta | str | None = None,
        origin: str | DatetimeLike = "start_day",
        restore_coord_dims: bool | None = None,
        **indexer_kwargs: ResampleCompatible | Resampler,
    ) -> DataArrayResample:
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

        Examples
        --------
        Downsample monthly time-series data to seasonal data:

        >>> da = xr.DataArray(
        ...     np.linspace(0, 11, num=12),
        ...     coords=[
        ...         pd.date_range(
        ...             "1999-12-15",
        ...             periods=12,
        ...             freq=pd.DateOffset(months=1),
        ...         )
        ...     ],
        ...     dims="time",
        ... )
        >>> da
        <xarray.DataArray (time: 12)> Size: 96B
        array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11.])
        Coordinates:
          * time     (time) datetime64[ns] 96B 1999-12-15 2000-01-15 ... 2000-11-15
        >>> da.resample(time="QS-DEC").mean()
        <xarray.DataArray (time: 4)> Size: 32B
        array([ 1.,  4.,  7., 10.])
        Coordinates:
          * time     (time) datetime64[ns] 32B 1999-12-01 2000-03-01 ... 2000-09-01

        Upsample monthly time-series data to daily data:

        >>> da.resample(time="1D").interpolate("linear")  # +doctest: ELLIPSIS
        <xarray.DataArray (time: 337)> Size: 3kB
        array([ 0.        ,  0.03225806,  0.06451613,  0.09677419,  0.12903226,
                0.16129032,  0.19354839,  0.22580645,  0.25806452,  0.29032258,
                0.32258065,  0.35483871,  0.38709677,  0.41935484,  0.4516129 ,
                0.48387097,  0.51612903,  0.5483871 ,  0.58064516,  0.61290323,
                0.64516129,  0.67741935,  0.70967742,  0.74193548,  0.77419355,
                0.80645161,  0.83870968,  0.87096774,  0.90322581,  0.93548387,
                0.96774194,  1.        ,  ...,
                9.        ,  9.03333333,  9.06666667,  9.1       ,  9.13333333,
                9.16666667,  9.2       ,  9.23333333,  9.26666667,  9.3       ,
                9.33333333,  9.36666667,  9.4       ,  9.43333333,  9.46666667,
                9.5       ,  9.53333333,  9.56666667,  9.6       ,  9.63333333,
                9.66666667,  9.7       ,  9.73333333,  9.76666667,  9.8       ,
                9.83333333,  9.86666667,  9.9       ,  9.93333333,  9.96666667,
               10.        , 10.03225806, 10.06451613, 10.09677419, 10.12903226,
               10.16129032, 10.19354839, 10.22580645, 10.25806452, 10.29032258,
               10.32258065, 10.35483871, 10.38709677, 10.41935484, 10.4516129 ,
               10.48387097, 10.51612903, 10.5483871 , 10.58064516, 10.61290323,
               10.64516129, 10.67741935, 10.70967742, 10.74193548, 10.77419355,
               10.80645161, 10.83870968, 10.87096774, 10.90322581, 10.93548387,
               10.96774194, 11.        ])
        Coordinates:
          * time     (time) datetime64[ns] 3kB 1999-12-15 1999-12-16 ... 2000-11-15

        Limit scope of upsampling method

        >>> da.resample(time="1D").nearest(tolerance="1D")
        <xarray.DataArray (time: 337)> Size: 3kB
        array([ 0.,  0., nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan,  1.,  1.,  1., nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan, nan, nan, nan, nan, nan,  2.,  2.,  2., nan,
               nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,  3.,
                3.,  3., nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan,  4.,  4.,  4., nan, nan, nan, nan, nan, ...,
               nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan, nan, 10., 10., 10., nan, nan, nan, nan, nan,
               nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,
               nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, 11., 11.])
        Coordinates:
          * time     (time) datetime64[ns] 3kB 1999-12-15 1999-12-16 ... 2000-11-15

        See Also
        --------
        Dataset.resample
        pandas.Series.resample
        pandas.DataFrame.resample

        References
        ----------
        .. [1] https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases
        """
        from xarray.core.resample import DataArrayResample

        return self._resample(
            resample_cls=DataArrayResample,
            indexer=indexer,
            skipna=skipna,
            closed=closed,
            label=label,
            offset=offset,
            origin=origin,
            restore_coord_dims=restore_coord_dims,
            **indexer_kwargs,
        )

    def to_dask_dataframe(
        self,
        dim_order: Sequence[Hashable] | None = None,
        set_index: bool = False,
    ) -> DaskDataFrame:
        """Convert this array into a dask.dataframe.DataFrame.

        Parameters
        ----------
        dim_order : Sequence of Hashable or None , optional
            Hierarchical dimension order for the resulting dataframe.
            Array content is transposed to this order and then written out as flat
            vectors in contiguous order, so the last dimension in this list
            will be contiguous in the resulting DataFrame. This has a major influence
            on which operations are efficient on the resulting dask dataframe.
        set_index : bool, default: False
            If set_index=True, the dask DataFrame is indexed by this dataset's
            coordinate. Since dask DataFrames do not support multi-indexes,
            set_index only works if the dataset only contains one dimension.

        Returns
        -------
        dask.dataframe.DataFrame

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.arange(4 * 2 * 2).reshape(4, 2, 2),
        ...     dims=("time", "lat", "lon"),
        ...     coords={
        ...         "time": np.arange(4),
        ...         "lat": [-30, -20],
        ...         "lon": [120, 130],
        ...     },
        ...     name="eg_dataarray",
        ...     attrs={"units": "Celsius", "description": "Random temperature data"},
        ... )
        >>> da.to_dask_dataframe(["lat", "lon", "time"]).compute()
            lat  lon  time  eg_dataarray
        0   -30  120     0             0
        1   -30  120     1             4
        2   -30  120     2             8
        3   -30  120     3            12
        4   -30  130     0             1
        5   -30  130     1             5
        6   -30  130     2             9
        7   -30  130     3            13
        8   -20  120     0             2
        9   -20  120     1             6
        10  -20  120     2            10
        11  -20  120     3            14
        12  -20  130     0             3
        13  -20  130     1             7
        14  -20  130     2            11
        15  -20  130     3            15
        """
        if self.name is None:
            raise ValueError(
                "Cannot convert an unnamed DataArray to a "
                "dask dataframe : use the ``.rename`` method to assign a name."
            )
        name = self.name
        ds = self._to_dataset_whole(name, shallow_copy=False)
        return ds.to_dask_dataframe(dim_order, set_index)

    # this needs to be at the end, or mypy will confuse with `str`
    # https://mypy.readthedocs.io/en/latest/common_issues.html#dealing-with-conflicting-names
    str = utils.UncachedAccessor(StringAccessor["DataArray"])

    def drop_attrs(self, *, deep: bool = True) -> Self:
        """
        Removes all attributes from the DataArray.

        Parameters
        ----------
        deep : bool, default True
            Removes attributes from coordinates.

        Returns
        -------
        DataArray
        """
        if not deep:
            return self._replace(attrs={})
        else:
            return (
                self._to_temp_dataset()
                .drop_attrs(deep=deep)
                .pipe(self._from_temp_dataset)
            )
