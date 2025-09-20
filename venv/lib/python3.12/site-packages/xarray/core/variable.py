from __future__ import annotations

import copy
import itertools
import math
import numbers
import warnings
from collections.abc import Callable, Hashable, Mapping, Sequence
from functools import partial
from types import EllipsisType
from typing import TYPE_CHECKING, Any, NoReturn, cast

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from packaging.version import Version

import xarray as xr  # only for Dataset and DataArray
from xarray.compat.array_api_compat import to_like_array
from xarray.computation import ops
from xarray.computation.arithmetic import VariableArithmetic
from xarray.core import common, dtypes, duck_array_ops, indexing, nputils, utils
from xarray.core.common import AbstractArray
from xarray.core.extension_array import PandasExtensionArray
from xarray.core.indexing import (
    BasicIndexer,
    CoordinateTransformIndexingAdapter,
    OuterIndexer,
    PandasIndexingAdapter,
    VectorizedIndexer,
    as_indexable,
)
from xarray.core.options import OPTIONS, _get_keep_attrs
from xarray.core.utils import (
    OrderedSet,
    _default,
    consolidate_dask_from_array_kwargs,
    decode_numpy_dict_values,
    drop_dims_from_indexers,
    either_dict_or_kwargs,
    emit_user_level_warning,
    ensure_us_time_resolution,
    infix_dims,
    is_allowed_extension_array,
    is_dict_like,
    is_duck_array,
    is_duck_dask_array,
    maybe_coerce_to_str,
)
from xarray.namedarray.core import NamedArray, _raise_if_any_duplicate_dimensions
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import (
    async_to_duck_array,
    integer_types,
    is_0d_dask_array,
    is_chunked_array,
    to_duck_array,
)
from xarray.namedarray.utils import module_available
from xarray.util.deprecation_helpers import _deprecate_positional_args, deprecate_dims

NON_NUMPY_SUPPORTED_ARRAY_TYPES = (
    indexing.ExplicitlyIndexed,
    pd.Index,
    pd.api.extensions.ExtensionArray,
    PandasExtensionArray,
)
# https://github.com/python/mypy/issues/224
BASIC_INDEXING_TYPES = integer_types + (slice,)
UNSUPPORTED_EXTENSION_ARRAY_TYPES = (
    pd.arrays.DatetimeArray,
    pd.arrays.TimedeltaArray,
    pd.arrays.NumpyExtensionArray,  # type: ignore[attr-defined]
)

if TYPE_CHECKING:
    from xarray.core.types import (
        Dims,
        ErrorOptionsWithWarn,
        PadModeOptions,
        PadReflectOptions,
        QuantileMethods,
        Self,
        T_Chunks,
        T_DuckArray,
        T_VarPadConstantValues,
    )
    from xarray.namedarray.parallelcompat import ChunkManagerEntrypoint


class MissingDimensionsError(ValueError):
    """Error class used when we can't safely guess a dimension name."""

    # inherits from ValueError for backward compatibility
    # TODO: move this to an xarray.exceptions module?


def as_variable(
    obj: T_DuckArray | Any, name=None, auto_convert: bool = True
) -> Variable | IndexVariable:
    """Convert an object into a Variable.

    Parameters
    ----------
    obj : object
        Object to convert into a Variable.

        - If the object is already a Variable, return a shallow copy.
        - Otherwise, if the object has 'dims' and 'data' attributes, convert
          it into a new Variable.
        - If all else fails, attempt to convert the object into a Variable by
          unpacking it into the arguments for creating a new Variable.
    name : str, optional
        If provided:

        - `obj` can be a 1D array, which is assumed to label coordinate values
          along a dimension of this given name.
        - Variables with name matching one of their dimensions are converted
          into `IndexVariable` objects.
    auto_convert : bool, optional
        For internal use only! If True, convert a "dimension" variable into
        an IndexVariable object (deprecated).

    Returns
    -------
    var : Variable
        The newly created variable.

    """
    from xarray.core.dataarray import DataArray

    # TODO: consider extending this method to automatically handle Iris and
    if isinstance(obj, DataArray):
        # extract the primary Variable from DataArrays
        obj = obj.variable

    if isinstance(obj, Variable):
        obj = obj.copy(deep=False)
    elif isinstance(obj, tuple):
        try:
            dims_, data_, *attrs = obj
        except ValueError as err:
            raise ValueError(
                f"Tuple {obj} is not in the form (dims, data[, attrs])"
            ) from err

        if isinstance(data_, DataArray):
            raise TypeError(
                f"Variable {name!r}: Using a DataArray object to construct a variable is"
                " ambiguous, please extract the data using the .data property."
            )
        try:
            obj = Variable(dims_, data_, *attrs)
        except (TypeError, ValueError) as error:
            raise error.__class__(
                f"Variable {name!r}: Could not convert tuple of form "
                f"(dims, data[, attrs, encoding]): {obj} to Variable."
            ) from error
    elif utils.is_scalar(obj):
        obj = Variable([], obj)
    elif isinstance(obj, pd.Index | IndexVariable) and obj.name is not None:
        obj = Variable(obj.name, obj)
    elif isinstance(obj, set | dict):
        raise TypeError(f"variable {name!r} has invalid type {type(obj)!r}")
    elif name is not None:
        data: T_DuckArray = as_compatible_data(obj)
        if data.ndim != 1:
            raise MissingDimensionsError(
                f"cannot set variable {name!r} with {data.ndim!r}-dimensional data "
                "without explicit dimension names. Pass a tuple of "
                "(dims, data) instead."
            )
        obj = Variable(name, data, fastpath=True)
    else:
        raise TypeError(
            f"Variable {name!r}: unable to convert object into a variable without an "
            f"explicit list of dimensions: {obj!r}"
        )

    if auto_convert and name is not None and name in obj.dims and obj.ndim == 1:
        # automatically convert the Variable into an Index
        emit_user_level_warning(
            f"variable {name!r} with name matching its dimension will not be "
            "automatically converted into an `IndexVariable` object in the future.",
            FutureWarning,
        )
        obj = obj.to_index_variable()

    return obj


def _maybe_wrap_data(data):
    """
    Put pandas.Index and numpy.ndarray arguments in adapter objects to ensure
    they can be indexed properly.

    NumpyArrayAdapter, PandasIndexingAdapter and LazilyIndexedArray should
    all pass through unmodified.
    """
    if isinstance(data, pd.Index):
        return PandasIndexingAdapter(data)
    if isinstance(data, UNSUPPORTED_EXTENSION_ARRAY_TYPES):
        return data.to_numpy()
    if isinstance(
        data, pd.api.extensions.ExtensionArray
    ) and is_allowed_extension_array(data):
        return PandasExtensionArray(data)
    return data


def _possibly_convert_objects(values):
    """Convert object arrays into datetime64 and timedelta64 according
    to the pandas convention.  For backwards compat, as of 3.0.0 pandas,
    object dtype inputs are cast to strings by `pandas.Series`
    but we output them as object dtype with the input metadata preserved as well.


    * datetime.datetime
    * datetime.timedelta
    * pd.Timestamp
    * pd.Timedelta
    """
    as_series = pd.Series(values.ravel(), copy=False)
    result = np.asarray(as_series).reshape(values.shape)
    if not result.flags.writeable:
        # GH8843, pandas copy-on-write mode creates read-only arrays by default
        try:
            result.flags.writeable = True
        except ValueError:
            result = result.copy()
    # For why we need this behavior: https://github.com/pandas-dev/pandas/issues/61938
    # Object datatype inputs that are strings
    # will be converted to strings by `pandas.Series`, and as of 3.0.0, lose
    # `dtype.metadata`.  If the roundtrip back to numpy in this function yields an
    # object array again, the dtype.metadata will be preserved.
    if (
        result.dtype.kind == "O"
        and values.dtype.kind == "O"
        and Version(pd.__version__) >= Version("3.0.0dev0")
    ):
        result.dtype = values.dtype
    return result


def as_compatible_data(
    data: T_DuckArray | ArrayLike, fastpath: bool = False
) -> T_DuckArray:
    """Prepare and wrap data to put in a Variable.

    - If data does not have the necessary attributes, convert it to ndarray.
    - If it's a pandas.Timestamp, convert it to datetime64.
    - If data is already a pandas or xarray object (other than an Index), just
      use the values.

    Finally, wrap it up with an adapter if necessary.
    """
    if fastpath and getattr(data, "ndim", None) is not None:
        return cast("T_DuckArray", data)

    from xarray.core.dataarray import DataArray

    # TODO: do this uwrapping in the Variable/NamedArray constructor instead.
    if isinstance(data, Variable):
        return cast("T_DuckArray", data._data)

    # TODO: do this uwrapping in the DataArray constructor instead.
    if isinstance(data, DataArray):
        return cast("T_DuckArray", data._variable._data)

    def convert_non_numpy_type(data):
        return cast("T_DuckArray", _maybe_wrap_data(data))

    if isinstance(data, NON_NUMPY_SUPPORTED_ARRAY_TYPES):
        return convert_non_numpy_type(data)

    if isinstance(data, tuple):
        data = utils.to_0d_object_array(data)

    # we don't want nested self-described arrays
    if isinstance(data, pd.Series | pd.DataFrame):
        if (
            isinstance(data, pd.Series)
            and is_allowed_extension_array(data.array)
            # Some datetime types are not allowed as well as backing Variable types
            and not isinstance(data.array, UNSUPPORTED_EXTENSION_ARRAY_TYPES)
        ):
            pandas_data = data.array
        else:
            pandas_data = data.values  # type: ignore[assignment]
        if isinstance(pandas_data, NON_NUMPY_SUPPORTED_ARRAY_TYPES):
            return convert_non_numpy_type(pandas_data)
        else:
            data = pandas_data

    if isinstance(data, np.ma.MaskedArray):
        mask = np.ma.getmaskarray(data)
        if mask.any():
            dtype, fill_value = dtypes.maybe_promote(data.dtype)
            data = duck_array_ops.where_method(data, ~mask, fill_value)
        else:
            data = np.asarray(data)

    if isinstance(data, np.matrix):
        data = np.asarray(data)

    # immediately return array-like types except `numpy.ndarray` and `numpy` scalars
    # compare types with `is` instead of `isinstance` to allow `numpy.ndarray` subclasses
    is_numpy = type(data) is np.ndarray or isinstance(data, np.generic)
    if not is_numpy and (
        hasattr(data, "__array_function__") or hasattr(data, "__array_namespace__")
    ):
        return cast("T_DuckArray", data)

    # anything left will be converted to `numpy.ndarray`, including `numpy` scalars
    data = np.asarray(data)

    if data.dtype.kind in "OMm":
        data = _possibly_convert_objects(data)
    return _maybe_wrap_data(data)


def _as_array_or_item(data):
    """Return the given values as a numpy array, or as an individual item if
    it's a 0d datetime64 or timedelta64 array.

    Importantly, this function does not copy data if it is already an ndarray -
    otherwise, it will not be possible to update Variable values in place.

    This function mostly exists because 0-dimensional ndarrays with
    dtype=datetime64 are broken :(
    https://github.com/numpy/numpy/issues/4337
    https://github.com/numpy/numpy/issues/7619

    TODO: remove this (replace with np.asarray) once these issues are fixed
    """
    data = np.asarray(data)
    if data.ndim == 0:
        kind = data.dtype.kind
        if kind in "mM":
            unit, _ = np.datetime_data(data.dtype)
            if kind == "M":
                data = np.datetime64(data, unit)
            elif kind == "m":
                data = np.timedelta64(data, unit)
    return data


class Variable(NamedArray, AbstractArray, VariableArithmetic):
    """A netcdf-like variable consisting of dimensions, data and attributes
    which describe a single Array. A single Variable object is not fully
    described outside the context of its parent Dataset (if you want such a
    fully described object, use a DataArray instead).

    The main functional difference between Variables and numpy arrays is that
    numerical operations on Variables implement array broadcasting by dimension
    name. For example, adding an Variable with dimensions `('time',)` to
    another Variable with dimensions `('space',)` results in a new Variable
    with dimensions `('time', 'space')`. Furthermore, numpy reduce operations
    like ``mean`` or ``sum`` are overwritten to take a "dimension" argument
    instead of an "axis".

    Variables are light-weight objects used as the building block for datasets.
    They are more primitive objects, so operations with them provide marginally
    higher performance than using DataArrays. However, manipulating data in the
    form of a Dataset or DataArray should almost always be preferred, because
    they can use more complete metadata in context of coordinate labels.
    """

    __slots__ = ("_attrs", "_data", "_dims", "_encoding")

    def __init__(
        self,
        dims,
        data: T_DuckArray | ArrayLike,
        attrs=None,
        encoding=None,
        fastpath=False,
    ):
        """
        Parameters
        ----------
        dims : str or sequence of str
            Name(s) of the the data dimension(s). Must be either a string (only
            for 1D data) or a sequence of strings with length equal to the
            number of dimensions.
        data : array_like
            Data array which supports numpy-like data access.
        attrs : dict_like or None, optional
            Attributes to assign to the new variable. If None (default), an
            empty attribute dictionary is initialized.
            (see FAQ, :ref:`approach to metadata`)
        encoding : dict_like or None, optional
            Dictionary specifying how to encode this array's data into a
            serialized format like netCDF4. Currently used keys (for netCDF)
            include '_FillValue', 'scale_factor', 'add_offset' and 'dtype'.
            Well-behaved code to serialize a Variable should ignore
            unrecognized encoding items.
        """
        super().__init__(
            dims=dims, data=as_compatible_data(data, fastpath=fastpath), attrs=attrs
        )

        self._encoding: dict[Any, Any] | None = None
        if encoding is not None:
            self.encoding = encoding

    def _new(
        self,
        dims=_default,
        data=_default,
        attrs=_default,
    ):
        dims_ = copy.copy(self._dims) if dims is _default else dims

        if attrs is _default:
            attrs_ = None if self._attrs is None else self._attrs.copy()
        else:
            attrs_ = attrs

        if data is _default:
            return type(self)(dims_, copy.copy(self._data), attrs_)
        else:
            cls_ = type(self)
            return cls_(dims_, data, attrs_)

    @property
    def _in_memory(self) -> bool:
        if isinstance(
            self._data, PandasIndexingAdapter | CoordinateTransformIndexingAdapter
        ):
            return self._data._in_memory

        return isinstance(
            self._data,
            np.ndarray | np.number | PandasExtensionArray,
        ) or (
            isinstance(self._data, indexing.MemoryCachedArray)
            and isinstance(self._data.array, indexing.NumpyIndexingAdapter)
        )

    @property
    def data(self):
        """
        The Variable's data as an array. The underlying array type
        (e.g. dask, sparse, pint) is preserved.

        See Also
        --------
        Variable.to_numpy
        Variable.as_numpy
        Variable.values
        """
        if isinstance(self._data, PandasExtensionArray):
            duck_array = self._data.array
        elif isinstance(self._data, indexing.ExplicitlyIndexed):
            duck_array = self._data.get_duck_array()
        elif is_duck_array(self._data):
            duck_array = self._data
        else:
            duck_array = self.values
        if isinstance(duck_array, PandasExtensionArray):
            # even though PandasExtensionArray is a duck array,
            # we should not return the PandasExtensionArray wrapper,
            # and instead return the underlying data.
            return duck_array.array
        return duck_array

    @data.setter  # type: ignore[override,unused-ignore]
    def data(self, data: T_DuckArray | ArrayLike) -> None:
        data = as_compatible_data(data)
        self._check_shape(data)
        self._data = data

    def astype(
        self,
        dtype,
        *,
        order=None,
        casting=None,
        subok=None,
        copy=None,
        keep_attrs=True,
    ) -> Self:
        """
        Copy of the Variable object, with data cast to a specified type.

        Parameters
        ----------
        dtype : str or dtype
            Typecode or data-type to which the array is cast.
        order : {'C', 'F', 'A', 'K'}, optional
            Controls the memory layout order of the result. ‘C’ means C order,
            ‘F’ means Fortran order, ‘A’ means ‘F’ order if all the arrays are
            Fortran contiguous, ‘C’ order otherwise, and ‘K’ means as close to
            the order the array elements appear in memory as possible.
        casting : {'no', 'equiv', 'safe', 'same_kind', 'unsafe'}, optional
            Controls what kind of data casting may occur.

            * 'no' means the data types should not be cast at all.
            * 'equiv' means only byte-order changes are allowed.
            * 'safe' means only casts which can preserve values are allowed.
            * 'same_kind' means only safe casts or casts within a kind,
              like float64 to float32, are allowed.
            * 'unsafe' means any data conversions may be done.
        subok : bool, optional
            If True, then sub-classes will be passed-through, otherwise the
            returned array will be forced to be a base-class array.
        copy : bool, optional
            By default, astype always returns a newly allocated array. If this
            is set to False and the `dtype` requirement is satisfied, the input
            array is returned instead of a copy.
        keep_attrs : bool, optional
            By default, astype keeps attributes. Set to False to remove
            attributes in the returned object.

        Returns
        -------
        out : same as object
            New object with data cast to the specified type.

        Notes
        -----
        The ``order``, ``casting``, ``subok`` and ``copy`` arguments are only passed
        through to the ``astype`` method of the underlying array when a value
        different than ``None`` is supplied.
        Make sure to only supply these arguments if the underlying array class
        supports them.

        See Also
        --------
        numpy.ndarray.astype
        dask.array.Array.astype
        sparse.COO.astype
        """
        from xarray.computation.apply_ufunc import apply_ufunc

        kwargs = dict(order=order, casting=casting, subok=subok, copy=copy)
        kwargs = {k: v for k, v in kwargs.items() if v is not None}

        return apply_ufunc(
            duck_array_ops.astype,
            self,
            dtype,
            kwargs=kwargs,
            keep_attrs=keep_attrs,
            dask="allowed",
        )

    def _dask_finalize(self, results, array_func, *args, **kwargs):
        data = array_func(results, *args, **kwargs)
        return Variable(self._dims, data, attrs=self._attrs, encoding=self._encoding)

    @property
    def values(self) -> np.ndarray:
        """The variable's data as a numpy.ndarray"""
        return _as_array_or_item(self._data)

    @values.setter
    def values(self, values):
        self.data = values

    def to_base_variable(self) -> Variable:
        """Return this variable as a base xarray.Variable"""
        return Variable(
            self._dims, self._data, self._attrs, encoding=self._encoding, fastpath=True
        )

    to_variable = utils.alias(to_base_variable, "to_variable")

    def to_index_variable(self) -> IndexVariable:
        """Return this variable as an xarray.IndexVariable"""
        return IndexVariable(
            self._dims, self._data, self._attrs, encoding=self._encoding, fastpath=True
        )

    to_coord = utils.alias(to_index_variable, "to_coord")

    def _to_index(self) -> pd.Index:
        return self.to_index_variable()._to_index()

    def to_index(self) -> pd.Index:
        """Convert this variable to a pandas.Index"""
        return self.to_index_variable().to_index()

    def to_dict(
        self, data: bool | str = "list", encoding: bool = False
    ) -> dict[str, Any]:
        """Dictionary representation of variable."""
        item: dict[str, Any] = {
            "dims": self.dims,
            "attrs": decode_numpy_dict_values(self.attrs),
        }
        if data is not False:
            if data in [True, "list"]:
                item["data"] = ensure_us_time_resolution(self.to_numpy()).tolist()
            elif data == "array":
                item["data"] = ensure_us_time_resolution(self.data)
            else:
                msg = 'data argument must be bool, "list", or "array"'
                raise ValueError(msg)

        else:
            item.update({"dtype": str(self.dtype), "shape": self.shape})

        if encoding:
            item["encoding"] = dict(self.encoding)

        return item

    def _item_key_to_tuple(self, key):
        if is_dict_like(key):
            return tuple(key.get(dim, slice(None)) for dim in self.dims)
        else:
            return key

    def _broadcast_indexes(self, key):
        """Prepare an indexing key for an indexing operation.

        Parameters
        ----------
        key : int, slice, array-like, dict or tuple of integer, slice and array-like
            Any valid input for indexing.

        Returns
        -------
        dims : tuple
            Dimension of the resultant variable.
        indexers : IndexingTuple subclass
            Tuple of integer, array-like, or slices to use when indexing
            self._data. The type of this argument indicates the type of
            indexing to perform, either basic, outer or vectorized.
        new_order : Optional[Sequence[int]]
            Optional reordering to do on the result of indexing. If not None,
            the first len(new_order) indexing should be moved to these
            positions.
        """
        key = self._item_key_to_tuple(key)  # key is a tuple
        # key is a tuple of full size
        key = indexing.expanded_indexer(key, self.ndim)
        # Convert a scalar Variable to a 0d-array
        key = tuple(
            k.data if isinstance(k, Variable) and k.ndim == 0 else k for k in key
        )
        # Convert a 0d numpy arrays to an integer
        # dask 0d arrays are passed through
        key = tuple(
            k.item() if isinstance(k, np.ndarray) and k.ndim == 0 else k for k in key
        )

        if all(
            (isinstance(k, BASIC_INDEXING_TYPES) and not isinstance(k, bool))
            for k in key
        ):
            return self._broadcast_indexes_basic(key)

        self._validate_indexers(key)
        # Detect it can be mapped as an outer indexer
        # If all key is unlabeled, or
        # key can be mapped as an OuterIndexer.
        if all(not isinstance(k, Variable) for k in key):
            return self._broadcast_indexes_outer(key)

        # If all key is 1-dimensional and there are no duplicate labels,
        # key can be mapped as an OuterIndexer.
        dims = []
        for k, d in zip(key, self.dims, strict=True):
            if isinstance(k, Variable):
                if len(k.dims) > 1:
                    return self._broadcast_indexes_vectorized(key)
                dims.append(k.dims[0])
            elif not isinstance(k, integer_types):
                dims.append(d)
        if len(set(dims)) == len(dims):
            return self._broadcast_indexes_outer(key)

        return self._broadcast_indexes_vectorized(key)

    def _broadcast_indexes_basic(self, key):
        dims = tuple(
            dim
            for k, dim in zip(key, self.dims, strict=True)
            if not isinstance(k, integer_types)
        )
        return dims, BasicIndexer(key), None

    def _validate_indexers(self, key):
        """Make sanity checks"""
        for dim, k in zip(self.dims, key, strict=True):
            if not isinstance(k, BASIC_INDEXING_TYPES):
                if not isinstance(k, Variable):
                    if not is_duck_array(k):
                        k = np.asarray(k)
                    if k.ndim > 1:
                        raise IndexError(
                            "Unlabeled multi-dimensional array cannot be "
                            f"used for indexing: {k}"
                        )
                if k.dtype.kind == "b":
                    if self.shape[self.get_axis_num(dim)] != len(k):
                        raise IndexError(
                            f"Boolean array size {len(k):d} is used to index array "
                            f"with shape {self.shape}."
                        )
                    if k.ndim > 1:
                        raise IndexError(
                            f"{k.ndim}-dimensional boolean indexing is not supported. "
                        )
                    if is_duck_dask_array(k.data):
                        raise KeyError(
                            "Indexing with a boolean dask array is not allowed. "
                            "This will result in a dask array of unknown shape. "
                            "Such arrays are unsupported by Xarray."
                            "Please compute the indexer first using .compute()"
                        )
                    if getattr(k, "dims", (dim,)) != (dim,):
                        raise IndexError(
                            "Boolean indexer should be unlabeled or on the "
                            "same dimension to the indexed array. Indexer is "
                            f"on {k.dims} but the target dimension is {dim}."
                        )

    def _broadcast_indexes_outer(self, key):
        # drop dim if k is integer or if k is a 0d dask array
        dims = tuple(
            k.dims[0] if isinstance(k, Variable) else dim
            for k, dim in zip(key, self.dims, strict=True)
            if (not isinstance(k, integer_types) and not is_0d_dask_array(k))
        )

        new_key = []
        for k in key:
            if isinstance(k, Variable):
                k = k.data
            if not isinstance(k, BASIC_INDEXING_TYPES):
                if not is_duck_array(k):
                    k = np.asarray(k)
                if k.size == 0:
                    # Slice by empty list; numpy could not infer the dtype
                    k = k.astype(int)
                elif k.dtype.kind == "b":
                    (k,) = np.nonzero(k)
            new_key.append(k)

        return dims, OuterIndexer(tuple(new_key)), None

    def _broadcast_indexes_vectorized(self, key):
        variables = []
        out_dims_set = OrderedSet()
        for dim, value in zip(self.dims, key, strict=True):
            if isinstance(value, slice):
                out_dims_set.add(dim)
            else:
                variable = (
                    value
                    if isinstance(value, Variable)
                    else as_variable(value, name=dim, auto_convert=False)
                )
                if variable.dims == (dim,):
                    variable = variable.to_index_variable()
                if variable.dtype.kind == "b":  # boolean indexing case
                    (variable,) = variable._nonzero()

                variables.append(variable)
                out_dims_set.update(variable.dims)

        variable_dims = set()
        for variable in variables:
            variable_dims.update(variable.dims)

        slices = []
        for i, (dim, value) in enumerate(zip(self.dims, key, strict=True)):
            if isinstance(value, slice):
                if dim in variable_dims:
                    # We only convert slice objects to variables if they share
                    # a dimension with at least one other variable. Otherwise,
                    # we can equivalently leave them as slices aknd transpose
                    # the result. This is significantly faster/more efficient
                    # for most array backends.
                    values = np.arange(*value.indices(self.sizes[dim]))
                    variables.insert(i - len(slices), Variable((dim,), values))
                else:
                    slices.append((i, value))

        try:
            variables = _broadcast_compat_variables(*variables)
        except ValueError as err:
            raise IndexError(f"Dimensions of indexers mismatch: {key}") from err

        out_key = [variable.data for variable in variables]
        out_dims = tuple(out_dims_set)
        slice_positions = set()
        for i, value in slices:
            out_key.insert(i, value)
            new_position = out_dims.index(self.dims[i])
            slice_positions.add(new_position)

        if slice_positions:
            new_order = [i for i in range(len(out_dims)) if i not in slice_positions]
        else:
            new_order = None

        return out_dims, VectorizedIndexer(tuple(out_key)), new_order

    def __getitem__(self, key) -> Self:
        """Return a new Variable object whose contents are consistent with
        getting the provided key from the underlying data.

        NB. __getitem__ and __setitem__ implement xarray-style indexing,
        where if keys are unlabeled arrays, we index the array orthogonally
        with them. If keys are labeled array (such as Variables), they are
        broadcasted with our usual scheme and then the array is indexed with
        the broadcasted key, like numpy's fancy indexing.

        If you really want to do indexing like `x[x > 0]`, manipulate the numpy
        array `x.values` directly.
        """
        dims, indexer, new_order = self._broadcast_indexes(key)
        indexable = as_indexable(self._data)

        data = indexing.apply_indexer(indexable, indexer)

        if new_order:
            data = duck_array_ops.moveaxis(data, range(len(new_order)), new_order)
        return self._finalize_indexing_result(dims, data)

    def _finalize_indexing_result(self, dims, data) -> Self:
        """Used by IndexVariable to return IndexVariable objects when possible."""
        return self._replace(dims=dims, data=data)

    def _getitem_with_mask(self, key, fill_value=dtypes.NA):
        """Index this Variable with -1 remapped to fill_value."""
        # TODO(shoyer): expose this method in public API somewhere (isel?) and
        # use it for reindex.
        # TODO(shoyer): add a sanity check that all other integers are
        # non-negative
        # TODO(shoyer): add an optimization, remapping -1 to an adjacent value
        # that is actually indexed rather than mapping it to the last value
        # along each axis.

        if fill_value is dtypes.NA:
            fill_value = dtypes.get_fill_value(self.dtype)

        dims, indexer, new_order = self._broadcast_indexes(key)

        if self.size:
            if is_duck_dask_array(self._data):
                # dask's indexing is faster this way; also vindex does not
                # support negative indices yet:
                # https://github.com/dask/dask/pull/2967
                actual_indexer = indexing.posify_mask_indexer(indexer)
            else:
                actual_indexer = indexer

            indexable = as_indexable(self._data)
            data = indexing.apply_indexer(indexable, actual_indexer)

            mask = indexing.create_mask(indexer, self.shape, data)
            # we need to invert the mask in order to pass data first. This helps
            # pint to choose the correct unit
            # TODO: revert after https://github.com/hgrecco/pint/issues/1019 is fixed
            mask = to_like_array(mask, data)
            data = duck_array_ops.where(
                duck_array_ops.logical_not(mask), data, fill_value
            )
        else:
            # array cannot be indexed along dimensions of size 0, so just
            # build the mask directly instead.
            mask = indexing.create_mask(indexer, self.shape)
            data = duck_array_ops.broadcast_to(fill_value, getattr(mask, "shape", ()))

        if new_order:
            data = duck_array_ops.moveaxis(data, range(len(new_order)), new_order)
        return self._finalize_indexing_result(dims, data)

    def __setitem__(self, key, value):
        """__setitem__ is overloaded to access the underlying numpy values with
        orthogonal indexing.

        See __getitem__ for more details.
        """
        dims, index_tuple, new_order = self._broadcast_indexes(key)

        if not isinstance(value, Variable):
            value = as_compatible_data(value)
            if value.ndim > len(dims):
                raise ValueError(
                    f"shape mismatch: value array of shape {value.shape} could not be "
                    f"broadcast to indexing result with {len(dims)} dimensions"
                )
            if value.ndim == 0:
                value = Variable((), value)
            else:
                value = Variable(dims[-value.ndim :], value)
        # broadcast to become assignable
        value = value.set_dims(dims).data

        if new_order:
            value = duck_array_ops.asarray(value)
            value = value[(len(dims) - value.ndim) * (np.newaxis,) + (Ellipsis,)]
            value = duck_array_ops.moveaxis(value, new_order, range(len(new_order)))

        indexable = as_indexable(self._data)
        indexing.set_with_indexer(indexable, index_tuple, value)

    @property
    def encoding(self) -> dict[Any, Any]:
        """Dictionary of encodings on this variable."""
        if self._encoding is None:
            encoding: dict[Any, Any] = {}
            self._encoding = encoding
        return self._encoding

    @encoding.setter
    def encoding(self, value):
        try:
            self._encoding = dict(value)
        except ValueError as err:
            raise ValueError("encoding must be castable to a dictionary") from err

    def reset_encoding(self) -> Self:
        warnings.warn(
            "reset_encoding is deprecated since 2023.11, use `drop_encoding` instead",
            stacklevel=2,
        )
        return self.drop_encoding()

    def drop_encoding(self) -> Self:
        """Return a new Variable without encoding."""
        return self._replace(encoding={})

    def _copy(
        self,
        deep: bool = True,
        data: T_DuckArray | ArrayLike | None = None,
        memo: dict[int, Any] | None = None,
    ) -> Self:
        if data is None:
            data_old = self._data

            if not isinstance(data_old, indexing.MemoryCachedArray):
                ndata = data_old
            else:
                # don't share caching between copies
                # TODO: MemoryCachedArray doesn't match the array api:
                ndata = indexing.MemoryCachedArray(data_old.array)  # type: ignore[assignment]

            if deep:
                ndata = copy.deepcopy(ndata, memo)

        else:
            ndata = as_compatible_data(data)
            if self.shape != ndata.shape:  # type: ignore[attr-defined]
                raise ValueError(
                    f"Data shape {ndata.shape} must match shape of object {self.shape}"  # type: ignore[attr-defined]
                )

        attrs = copy.deepcopy(self._attrs, memo) if deep else copy.copy(self._attrs)
        encoding = (
            copy.deepcopy(self._encoding, memo) if deep else copy.copy(self._encoding)
        )

        # note: dims is already an immutable tuple
        return self._replace(data=ndata, attrs=attrs, encoding=encoding)

    def _replace(
        self,
        dims=_default,
        data=_default,
        attrs=_default,
        encoding=_default,
    ) -> Self:
        if dims is _default:
            dims = copy.copy(self._dims)
        if data is _default:
            data = copy.copy(self.data)
        if attrs is _default:
            attrs = copy.copy(self._attrs)
        if encoding is _default:
            encoding = copy.copy(self._encoding)
        return type(self)(dims, data, attrs, encoding, fastpath=True)

    def load(self, **kwargs) -> Self:
        """Trigger loading data into memory and return this variable.

        Data will be computed and/or loaded from disk or a remote source.

        Unlike ``.compute``, the original variable is modified and returned.

        Normally, it should not be necessary to call this method in user code,
        because all xarray functions should either work on deferred data or
        load data automatically.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.array.compute``.

        Returns
        -------
        object : Variable
            Same object but with lazy data as an in-memory array.

        See Also
        --------
        dask.array.compute
        Variable.compute
        Variable.load_async
        DataArray.load
        Dataset.load
        """
        self._data = to_duck_array(self._data, **kwargs)
        return self

    async def load_async(self, **kwargs) -> Self:
        """Trigger and await asynchronous loading of data into memory and return this variable.

        Data will be computed and/or loaded from disk or a remote source.

        Unlike ``.compute``, the original variable is modified and returned.

        Only works when opening data lazily from IO storage backends which support lazy asynchronous loading.
        Otherwise will raise a NotImplementedError.

        Note users are expected to limit concurrency themselves - xarray does not internally limit concurrency in any way.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.array.compute``.

        Returns
        -------
        object : Variable
            Same object but with lazy data as an in-memory array.

        See Also
        --------
        dask.array.compute
        Variable.load
        Variable.compute
        DataArray.load_async
        Dataset.load_async
        """
        self._data = await async_to_duck_array(self._data, **kwargs)
        return self

    def compute(self, **kwargs) -> Self:
        """Trigger loading data into memory and return a new variable.

        Data will be computed and/or loaded from disk or a remote source.

        The original variable is left unaltered.

        Normally, it should not be necessary to call this method in user code,
        because all xarray functions should either work on deferred data or
        load data automatically.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments passed on to ``dask.array.compute``.

        Returns
        -------
        object : Variable
            New object with the data as an in-memory array.

        See Also
        --------
        dask.array.compute
        Variable.load
        Variable.load_async
        DataArray.compute
        Dataset.compute
        """
        new = self.copy(deep=False)
        return new.load(**kwargs)

    def _shuffle(
        self, indices: list[list[int]], dim: Hashable, chunks: T_Chunks
    ) -> Self:
        # TODO (dcherian): consider making this public API
        array = self._data
        if is_chunked_array(array):
            chunkmanager = get_chunked_array_type(array)
            return self._replace(
                data=chunkmanager.shuffle(
                    array,
                    indexer=indices,
                    axis=self.get_axis_num(dim),
                    chunks=chunks,
                )
            )
        else:
            return self.isel({dim: np.concatenate(indices)})

    def isel(
        self,
        indexers: Mapping[Any, Any] | None = None,
        missing_dims: ErrorOptionsWithWarn = "raise",
        **indexers_kwargs: Any,
    ) -> Self:
        """Return a new array indexed along the specified dimension(s).

        Parameters
        ----------
        **indexers : {dim: indexer, ...}
            Keyword arguments with names matching dimensions and values given
            by integers, slice objects or arrays.
        missing_dims : {"raise", "warn", "ignore"}, default: "raise"
            What to do if dimensions that should be selected from are not present in the
            DataArray:
            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        Returns
        -------
        obj : Array object
            A new Array with the selected data and dimensions. In general,
            the new variable's data will be a view of this variable's data,
            unless numpy fancy indexing was triggered by using an array
            indexer, in which case the data will be a copy.
        """
        indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "isel")

        indexers = drop_dims_from_indexers(indexers, self.dims, missing_dims)

        key = tuple(indexers.get(dim, slice(None)) for dim in self.dims)
        return self[key]

    def squeeze(self, dim=None):
        """Return a new object with squeezed data.

        Parameters
        ----------
        dim : None or str or tuple of str, optional
            Selects a subset of the length one dimensions. If a dimension is
            selected with length greater than one, an error is raised. If
            None, all length one dimensions are squeezed.

        Returns
        -------
        squeezed : same type as caller
            This object, but with with all or a subset of the dimensions of
            length 1 removed.

        See Also
        --------
        numpy.squeeze
        """
        dims = common.get_squeeze_dims(self, dim)
        return self.isel(dict.fromkeys(dims, 0))

    def _shift_one_dim(self, dim, count, fill_value=dtypes.NA):
        axis = self.get_axis_num(dim)

        if count > 0:
            keep = slice(None, -count)
        elif count < 0:
            keep = slice(-count, None)
        else:
            keep = slice(None)

        trimmed_data = self[(slice(None),) * axis + (keep,)].data

        if fill_value is dtypes.NA:
            dtype, fill_value = dtypes.maybe_promote(self.dtype)
        else:
            dtype = self.dtype

        width = min(abs(count), self.shape[axis])
        dim_pad = (width, 0) if count >= 0 else (0, width)
        pads = [(0, 0) if d != dim else dim_pad for d in self.dims]

        data = duck_array_ops.pad(
            duck_array_ops.astype(trimmed_data, dtype),
            pads,
            mode="constant",
            constant_values=fill_value,
        )

        if is_duck_dask_array(data):
            # chunked data should come out with the same chunks; this makes
            # it feasible to combine shifted and unshifted data
            # TODO: remove this once dask.array automatically aligns chunks
            data = data.rechunk(self.data.chunks)

        return self._replace(data=data)

    def shift(self, shifts=None, fill_value=dtypes.NA, **shifts_kwargs):
        """
        Return a new Variable with shifted data.

        Parameters
        ----------
        shifts : mapping of the form {dim: offset}
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
        shifted : Variable
            Variable with the same dimensions and attributes but shifted data.
        """
        shifts = either_dict_or_kwargs(shifts, shifts_kwargs, "shift")
        result = self
        for dim, count in shifts.items():
            result = result._shift_one_dim(dim, count, fill_value=fill_value)
        return result

    def _pad_options_dim_to_index(
        self,
        pad_option: Mapping[Any, int | float | tuple[int, int] | tuple[float, float]],
        fill_with_shape=False,
    ):
        # change number values to a tuple of two of those values
        for k, v in pad_option.items():
            if isinstance(v, numbers.Number):
                pad_option[k] = (v, v)

        if fill_with_shape:
            return [
                pad_option.get(d, (n, n))
                for d, n in zip(self.dims, self.data.shape, strict=True)
            ]
        return [pad_option.get(d, (0, 0)) for d in self.dims]

    def pad(
        self,
        pad_width: Mapping[Any, int | tuple[int, int]] | None = None,
        mode: PadModeOptions = "constant",
        stat_length: (
            int | tuple[int, int] | Mapping[Any, tuple[int, int]] | None
        ) = None,
        constant_values: T_VarPadConstantValues | None = None,
        end_values: int | tuple[int, int] | Mapping[Any, tuple[int, int]] | None = None,
        reflect_type: PadReflectOptions = None,
        keep_attrs: bool | None = None,
        **pad_width_kwargs: Any,
    ):
        """
        Return a new Variable with padded data.

        Parameters
        ----------
        pad_width : mapping of hashable to tuple of int
            Mapping with the form of {dim: (pad_before, pad_after)}
            describing the number of values padded along each dimension.
            {dim: pad} is a shortcut for pad_before = pad_after = pad
        mode : str, default: "constant"
            See numpy / Dask docs
        stat_length : int, tuple or mapping of hashable to tuple
            Used in 'maximum', 'mean', 'median', and 'minimum'.  Number of
            values at edge of each axis used to calculate the statistic value.
        constant_values : scalar, tuple or mapping of hashable to scalar or tuple
            Used in 'constant'.  The values to set the padded values for each
            axis.
        end_values : scalar, tuple or mapping of hashable to tuple
            Used in 'linear_ramp'.  The values used for the ending value of the
            linear_ramp and that will form the edge of the padded array.
        reflect_type : {"even", "odd"}, optional
            Used in "reflect", and "symmetric".  The "even" style is the
            default with an unaltered reflection around the edge value.  For
            the "odd" style, the extended part of the array is created by
            subtracting the reflected values from two times the edge value.
        keep_attrs : bool, optional
            If True, the variable's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        **pad_width_kwargs
            One of pad_width or pad_width_kwargs must be provided.

        Returns
        -------
        padded : Variable
            Variable with the same dimensions and attributes but padded data.
        """
        pad_width = either_dict_or_kwargs(pad_width, pad_width_kwargs, "pad")

        # change default behaviour of pad with mode constant
        if mode == "constant" and (
            constant_values is None or constant_values is dtypes.NA
        ):
            dtype, constant_values = dtypes.maybe_promote(self.dtype)
        else:
            dtype = self.dtype

        # create pad_options_kwargs, numpy requires only relevant kwargs to be nonempty
        if isinstance(stat_length, dict):
            stat_length = self._pad_options_dim_to_index(
                stat_length, fill_with_shape=True
            )
        if isinstance(constant_values, dict):
            constant_values = self._pad_options_dim_to_index(constant_values)
        if isinstance(end_values, dict):
            end_values = self._pad_options_dim_to_index(end_values)

        # workaround for bug in Dask's default value of stat_length https://github.com/dask/dask/issues/5303
        if stat_length is None and mode in ["maximum", "mean", "median", "minimum"]:
            stat_length = [(n, n) for n in self.data.shape]  # type: ignore[assignment]

        pad_width_by_index = self._pad_options_dim_to_index(pad_width)

        # create pad_options_kwargs, numpy/dask requires only relevant kwargs to be nonempty
        pad_option_kwargs: dict[str, Any] = {}
        if stat_length is not None:
            pad_option_kwargs["stat_length"] = stat_length
        if constant_values is not None:
            pad_option_kwargs["constant_values"] = constant_values
        if end_values is not None:
            pad_option_kwargs["end_values"] = end_values
        if reflect_type is not None:
            pad_option_kwargs["reflect_type"] = reflect_type

        array = duck_array_ops.pad(
            duck_array_ops.astype(self.data, dtype, copy=False),
            pad_width_by_index,
            mode=mode,
            **pad_option_kwargs,
        )

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)
        attrs = self._attrs if keep_attrs else None

        return type(self)(self.dims, array, attrs=attrs)

    def _roll_one_dim(self, dim, count):
        axis = self.get_axis_num(dim)

        count %= self.shape[axis]
        if count != 0:
            indices = [slice(-count, None), slice(None, -count)]
        else:
            indices = [slice(None)]

        arrays = [self[(slice(None),) * axis + (idx,)].data for idx in indices]

        data = duck_array_ops.concatenate(arrays, axis)

        if is_duck_dask_array(data):
            # chunked data should come out with the same chunks; this makes
            # it feasible to combine shifted and unshifted data
            # TODO: remove this once dask.array automatically aligns chunks
            data = data.rechunk(self.data.chunks)

        return self._replace(data=data)

    def roll(self, shifts=None, **shifts_kwargs):
        """
        Return a new Variable with rolld data.

        Parameters
        ----------
        shifts : mapping of hashable to int
            Integer offset to roll along each of the given dimensions.
            Positive offsets roll to the right; negative offsets roll to the
            left.
        **shifts_kwargs
            The keyword arguments form of ``shifts``.
            One of shifts or shifts_kwargs must be provided.

        Returns
        -------
        shifted : Variable
            Variable with the same dimensions and attributes but rolled data.
        """
        shifts = either_dict_or_kwargs(shifts, shifts_kwargs, "roll")

        result = self
        for dim, count in shifts.items():
            result = result._roll_one_dim(dim, count)
        return result

    @deprecate_dims
    def transpose(
        self,
        *dim: Hashable | EllipsisType,
        missing_dims: ErrorOptionsWithWarn = "raise",
    ) -> Self:
        """Return a new Variable object with transposed dimensions.

        Parameters
        ----------
        *dim : Hashable, optional
            By default, reverse the dimensions. Otherwise, reorder the
            dimensions to this order.
        missing_dims : {"raise", "warn", "ignore"}, default: "raise"
            What to do if dimensions that should be selected from are not present in the
            Variable:
            - "raise": raise an exception
            - "warn": raise a warning, and ignore the missing dimensions
            - "ignore": ignore the missing dimensions

        Returns
        -------
        transposed : Variable
            The returned object has transposed data and dimensions with the
            same attributes as the original.

        Notes
        -----
        This operation returns a view of this variable's data. It is
        lazy for dask-backed Variables but not for numpy-backed Variables.

        See Also
        --------
        numpy.transpose
        """
        if len(dim) == 0:
            dim = self.dims[::-1]
        else:
            dim = tuple(infix_dims(dim, self.dims, missing_dims))

        if len(dim) < 2 or dim == self.dims:
            # no need to transpose if only one dimension
            # or dims are in same order
            return self.copy(deep=False)

        axes = self.get_axis_num(dim)
        data = as_indexable(self._data).transpose(axes)
        return self._replace(dims=dim, data=data)

    @property
    def T(self) -> Self:
        return self.transpose()

    @deprecate_dims
    def set_dims(self, dim, shape=None):
        """Return a new variable with given set of dimensions.
        This method might be used to attach new dimension(s) to variable.

        When possible, this operation does not copy this variable's data.

        Parameters
        ----------
        dim : str or sequence of str or dict
            Dimensions to include on the new variable. If a dict, values are
            used to provide the sizes of new dimensions; otherwise, new
            dimensions are inserted with length 1.

        Returns
        -------
        Variable
        """
        if isinstance(dim, str):
            dim = [dim]

        if shape is None and is_dict_like(dim):
            shape = tuple(dim.values())

        missing_dims = set(self.dims) - set(dim)
        if missing_dims:
            raise ValueError(
                f"new dimensions {dim!r} must be a superset of "
                f"existing dimensions {self.dims!r}"
            )

        self_dims = set(self.dims)
        expanded_dims = tuple(d for d in dim if d not in self_dims) + self.dims

        if self.dims == expanded_dims:
            # don't use broadcast_to unless necessary so the result remains
            # writeable if possible
            expanded_data = self.data
        elif shape is None or all(
            s == 1 for s, e in zip(shape, dim, strict=True) if e not in self_dims
        ):
            # "Trivial" broadcasting, i.e. simply inserting a new dimension
            # This is typically easier for duck arrays to implement
            # than the full "broadcast_to" semantics
            indexer = (None,) * (len(expanded_dims) - self.ndim) + (...,)
            expanded_data = self.data[indexer]
        else:  # elif shape is not None:
            dims_map = dict(zip(dim, shape, strict=True))
            tmp_shape = tuple(dims_map[d] for d in expanded_dims)
            expanded_data = duck_array_ops.broadcast_to(self._data, tmp_shape)

        expanded_var = Variable(
            expanded_dims, expanded_data, self._attrs, self._encoding, fastpath=True
        )
        return expanded_var.transpose(*dim)

    def _stack_once(self, dim: list[Hashable], new_dim: Hashable):
        if not set(dim) <= set(self.dims):
            raise ValueError(f"invalid existing dimensions: {dim}")

        if new_dim in self.dims:
            raise ValueError(
                "cannot create a new dimension with the same "
                "name as an existing dimension"
            )

        if len(dim) == 0:
            # don't stack
            return self.copy(deep=False)

        other_dims = [d for d in self.dims if d not in dim]
        dim_order = other_dims + list(dim)
        reordered = self.transpose(*dim_order)

        new_shape = reordered.shape[: len(other_dims)] + (-1,)
        new_data = duck_array_ops.reshape(reordered.data, new_shape)
        new_dims = reordered.dims[: len(other_dims)] + (new_dim,)

        return type(self)(
            new_dims, new_data, self._attrs, self._encoding, fastpath=True
        )

    @partial(deprecate_dims, old_name="dimensions")
    def stack(self, dim=None, **dim_kwargs):
        """
        Stack any number of existing dim into a single new dimension.

        New dim will be added at the end, and the order of the data
        along each new dimension will be in contiguous (C) order.

        Parameters
        ----------
        dim : mapping of hashable to tuple of hashable
            Mapping of form new_name=(dim1, dim2, ...) describing the
            names of new dim, and the existing dim that
            they replace.
        **dim_kwargs
            The keyword arguments form of ``dim``.
            One of dim or dim_kwargs must be provided.

        Returns
        -------
        stacked : Variable
            Variable with the same attributes but stacked data.

        See Also
        --------
        Variable.unstack
        """
        dim = either_dict_or_kwargs(dim, dim_kwargs, "stack")
        result = self
        for new_dim, dims in dim.items():
            result = result._stack_once(dims, new_dim)
        return result

    def _unstack_once_full(self, dim: Mapping[Any, int], old_dim: Hashable) -> Self:
        """
        Unstacks the variable without needing an index.

        Unlike `_unstack_once`, this function requires the existing dimension to
        contain the full product of the new dimensions.
        """
        new_dim_names = tuple(dim.keys())
        new_dim_sizes = tuple(dim.values())

        if old_dim not in self.dims:
            raise ValueError(f"invalid existing dimension: {old_dim}")

        if set(new_dim_names).intersection(self.dims):
            raise ValueError(
                "cannot create a new dimension with the same "
                "name as an existing dimension"
            )

        if math.prod(new_dim_sizes) != self.sizes[old_dim]:
            raise ValueError(
                "the product of the new dimension sizes must "
                "equal the size of the old dimension"
            )

        other_dims = [d for d in self.dims if d != old_dim]
        dim_order = other_dims + [old_dim]
        reordered = self.transpose(*dim_order)

        new_shape = reordered.shape[: len(other_dims)] + new_dim_sizes
        new_data = duck_array_ops.reshape(reordered.data, new_shape)
        new_dims = reordered.dims[: len(other_dims)] + new_dim_names

        return type(self)(
            new_dims, new_data, self._attrs, self._encoding, fastpath=True
        )

    def _unstack_once(
        self,
        index: pd.MultiIndex,
        dim: Hashable,
        fill_value=dtypes.NA,
        sparse: bool = False,
    ) -> Variable:
        """
        Unstacks this variable given an index to unstack and the name of the
        dimension to which the index refers.
        """

        reordered = self.transpose(..., dim)

        new_dim_sizes = [lev.size for lev in index.levels]
        new_dim_names = index.names
        indexer = index.codes

        # Potentially we could replace `len(other_dims)` with just `-1`
        other_dims = [d for d in self.dims if d != dim]
        new_shape = tuple(list(reordered.shape[: len(other_dims)]) + new_dim_sizes)
        new_dims = reordered.dims[: len(other_dims)] + tuple(new_dim_names)

        create_template: Callable
        if fill_value is dtypes.NA:
            is_missing_values = math.prod(new_shape) > math.prod(self.shape)
            if is_missing_values:
                dtype, fill_value = dtypes.maybe_promote(self.dtype)

                create_template = partial(
                    duck_array_ops.full_like, fill_value=fill_value
                )
            else:
                dtype = self.dtype
                fill_value = dtypes.get_fill_value(dtype)
                create_template = duck_array_ops.empty_like
        else:
            dtype = self.dtype
            create_template = partial(duck_array_ops.full_like, fill_value=fill_value)

        if sparse:
            # unstacking a dense multitindexed array to a sparse array
            from sparse import COO

            codes = zip(*index.codes, strict=True)
            if reordered.ndim == 1:
                indexes = codes
            else:
                sizes = itertools.product(*[range(s) for s in reordered.shape[:-1]])
                tuple_indexes = itertools.product(sizes, codes)
                indexes = (list(itertools.chain(*x)) for x in tuple_indexes)  # type: ignore[assignment]

            data = COO(
                coords=np.array(list(indexes)).T,
                data=self.data.astype(dtype).ravel(),
                fill_value=fill_value,
                shape=new_shape,
                sorted=index.is_monotonic_increasing,
            )

        else:
            data = create_template(self.data, shape=new_shape, dtype=dtype)

            # Indexer is a list of lists of locations. Each list is the locations
            # on the new dimension. This is robust to the data being sparse; in that
            # case the destinations will be NaN / zero.
            data[(..., *indexer)] = reordered

        return self.to_base_variable()._replace(dims=new_dims, data=data)

    @partial(deprecate_dims, old_name="dimensions")
    def unstack(self, dim=None, **dim_kwargs) -> Variable:
        """
        Unstack an existing dimension into multiple new dimensions.

        New dimensions will be added at the end, and the order of the data
        along each new dimension will be in contiguous (C) order.

        Note that unlike ``DataArray.unstack`` and ``Dataset.unstack``, this
        method requires the existing dimension to contain the full product of
        the new dimensions.

        Parameters
        ----------
        dim : mapping of hashable to mapping of hashable to int
            Mapping of the form old_dim={dim1: size1, ...} describing the
            names of existing dimensions, and the new dimensions and sizes
            that they map to.
        **dim_kwargs
            The keyword arguments form of ``dim``.
            One of dim or dim_kwargs must be provided.

        Returns
        -------
        unstacked : Variable
            Variable with the same attributes but unstacked data.

        See Also
        --------
        Variable.stack
        DataArray.unstack
        Dataset.unstack
        """
        dim = either_dict_or_kwargs(dim, dim_kwargs, "unstack")
        result = self
        for old_dim, dims in dim.items():
            result = result._unstack_once_full(dims, old_dim)
        return result

    def fillna(self, value):
        return ops.fillna(self, value)

    def where(self, cond, other=dtypes.NA):
        return ops.where_method(self, cond, other)

    def clip(self, min=None, max=None):
        """
        Return an array whose values are limited to ``[min, max]``.
        At least one of max or min must be given.

        Refer to `numpy.clip` for full documentation.

        See Also
        --------
        numpy.clip : equivalent function
        """
        from xarray.computation.apply_ufunc import apply_ufunc

        xp = duck_array_ops.get_array_namespace(self.data)
        return apply_ufunc(xp.clip, self, min, max, dask="allowed")

    def reduce(  # type: ignore[override]
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs,
    ) -> Variable:
        """Reduce this array by applying `func` along some dimension(s).

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `func(x, axis=axis, **kwargs)` to return the result of reducing an
            np.ndarray over an integer valued axis.
        dim : "...", str, Iterable of Hashable or None, optional
            Dimension(s) over which to apply `func`. By default `func` is
            applied over all dimensions.
        axis : int or Sequence of int, optional
            Axis(es) over which to apply `func`. Only one of the 'dim'
            and 'axis' arguments can be supplied. If neither are supplied, then
            the reduction is calculated over the flattened array (by calling
            `func(x)` without an axis argument).
        keep_attrs : bool, optional
            If True, the variable's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        keepdims : bool, default: False
            If True, the dimensions which are reduced are left in the result
            as dimensions of size one
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : Array
            Array with summarized data and the indicated dimension(s)
            removed.
        """
        keep_attrs_ = (
            _get_keep_attrs(default=False) if keep_attrs is None else keep_attrs
        )

        # Note that the call order for Variable.mean is
        #    Variable.mean -> NamedArray.mean -> Variable.reduce
        #    -> NamedArray.reduce
        result = super().reduce(
            func=func, dim=dim, axis=axis, keepdims=keepdims, **kwargs
        )

        # return Variable always to support IndexVariable
        return Variable(
            result.dims, result._data, attrs=result._attrs if keep_attrs_ else None
        )

    @classmethod
    def concat(
        cls,
        variables,
        dim="concat_dim",
        positions=None,
        shortcut=False,
        combine_attrs="override",
    ):
        """Concatenate variables along a new or existing dimension.

        Parameters
        ----------
        variables : iterable of Variable
            Arrays to stack together. Each variable is expected to have
            matching dimensions and shape except for along the stacked
            dimension.
        dim : str or DataArray, optional
            Name of the dimension to stack along. This can either be a new
            dimension name, in which case it is added along axis=0, or an
            existing dimension name, in which case the location of the
            dimension is unchanged. Where to insert the new dimension is
            determined by the first variable.
        positions : None or list of array-like, optional
            List of integer arrays which specifies the integer positions to
            which to assign each dataset along the concatenated dimension.
            If not supplied, objects are concatenated in the provided order.
        shortcut : bool, optional
            This option is used internally to speed-up groupby operations.
            If `shortcut` is True, some checks of internal consistency between
            arrays to concatenate are skipped.
        combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts", \
                         "override"}, default: "override"
            String indicating how to combine attrs of the objects being merged:

            - "drop": empty attrs on returned Dataset.
            - "identical": all attrs must be the same on every object.
            - "no_conflicts": attrs from all objects are combined, any that have
              the same name must also have the same value.
            - "drop_conflicts": attrs from all objects are combined, any that have
              the same name but different values are dropped.
            - "override": skip comparing and copy attrs from the first dataset to
              the result.

        Returns
        -------
        stacked : Variable
            Concatenated Variable formed by stacking all the supplied variables
            along the given dimension.
        """
        from xarray.structure.merge import merge_attrs

        if not isinstance(dim, str):
            (dim,) = dim.dims

        # can't do this lazily: we need to loop through variables at least
        # twice
        variables = list(variables)
        first_var = variables[0]
        first_var_dims = first_var.dims

        arrays = [v._data for v in variables]

        if dim in first_var_dims:
            axis = first_var.get_axis_num(dim)
            dims = first_var_dims
            data = duck_array_ops.concatenate(arrays, axis=axis)
            if positions is not None:
                # TODO: deprecate this option -- we don't need it for groupby
                # any more.
                indices = nputils.inverse_permutation(np.concatenate(positions))
                data = duck_array_ops.take(data, indices, axis=axis)
        else:
            axis = 0
            dims = (dim,) + first_var_dims
            data = duck_array_ops.stack(arrays, axis=axis)

        attrs = merge_attrs(
            [var.attrs for var in variables], combine_attrs=combine_attrs
        )
        encoding = dict(first_var.encoding)
        if not shortcut:
            for var in variables:
                if var.dims != first_var_dims:
                    raise ValueError(
                        f"Variable has dimensions {tuple(var.dims)} but first Variable has dimensions {tuple(first_var_dims)}"
                    )

        return cls(dims, data, attrs, encoding, fastpath=True)

    def equals(self, other, equiv=duck_array_ops.array_equiv):
        """True if two Variables have the same dimensions and values;
        otherwise False.

        Variables can still be equal (like pandas objects) if they have NaN
        values in the same locations.

        This method is necessary because `v1 == v2` for Variables
        does element-wise comparisons (like numpy.ndarrays).
        """
        other = getattr(other, "variable", other)
        try:
            return self.dims == other.dims and (
                self._data is other._data or equiv(self.data, other.data)
            )
        except (TypeError, AttributeError):
            return False

    def broadcast_equals(self, other, equiv=duck_array_ops.array_equiv):
        """True if two Variables have the values after being broadcast against
        each other; otherwise False.

        Variables can still be equal (like pandas objects) if they have NaN
        values in the same locations.
        """
        try:
            self, other = broadcast_variables(self, other)
        except (ValueError, AttributeError):
            return False
        return self.equals(other, equiv=equiv)

    def identical(self, other, equiv=duck_array_ops.array_equiv):
        """Like equals, but also checks attributes."""
        try:
            return utils.dict_equiv(self.attrs, other.attrs) and self.equals(
                other, equiv=equiv
            )
        except (TypeError, AttributeError):
            return False

    def no_conflicts(self, other, equiv=duck_array_ops.array_notnull_equiv):
        """True if the intersection of two Variable's non-null data is
        equal; otherwise false.

        Variables can thus still be equal if there are locations where either,
        or both, contain NaN values.
        """
        return self.broadcast_equals(other, equiv=equiv)

    def quantile(
        self,
        q: ArrayLike,
        dim: str | Sequence[Hashable] | None = None,
        method: QuantileMethods = "linear",
        keep_attrs: bool | None = None,
        skipna: bool | None = None,
        interpolation: QuantileMethods | None = None,
    ) -> Self:
        """Compute the qth quantile of the data along the specified dimension.

        Returns the qth quantiles(s) of the array elements.

        Parameters
        ----------
        q : float or sequence of float
            Quantile to compute, which must be between 0 and 1
            inclusive.
        dim : str or sequence of str, optional
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
            If True, the variable's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        skipna : bool, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or skipna=True has not been
            implemented (object, datetime64 or timedelta64).

        Returns
        -------
        quantiles : Variable
            If `q` is a single quantile, then the result
            is a scalar. If multiple percentiles are given, first axis of
            the result corresponds to the quantile and a quantile dimension
            is added to the return array. The other dimensions are the
            dimensions that remain after the reduction of the array.

        See Also
        --------
        numpy.nanquantile, pandas.Series.quantile, Dataset.quantile
        DataArray.quantile

        References
        ----------
        .. [1] R. J. Hyndman and Y. Fan,
           "Sample quantiles in statistical packages,"
           The American Statistician, 50(4), pp. 361-365, 1996
        """

        from xarray.computation.apply_ufunc import apply_ufunc

        if interpolation is not None:
            warnings.warn(
                "The `interpolation` argument to quantile was renamed to `method`.",
                FutureWarning,
                stacklevel=2,
            )

            if method != "linear":
                raise TypeError("Cannot pass interpolation and method keywords!")

            method = interpolation

        if skipna or (skipna is None and self.dtype.kind in "cfO"):
            _quantile_func = nputils.nanquantile
        else:
            _quantile_func = duck_array_ops.quantile

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)

        scalar = utils.is_scalar(q)
        q = np.atleast_1d(np.asarray(q, dtype=np.float64))

        if dim is None:
            dim = self.dims

        if utils.is_scalar(dim):
            dim = [dim]

        xp = duck_array_ops.get_array_namespace(self.data)

        def _wrapper(npa, **kwargs):
            # move quantile axis to end. required for apply_ufunc
            return xp.moveaxis(_quantile_func(npa, **kwargs), 0, -1)

        # jax requires hashable
        axis = tuple(range(-1, -1 * len(dim) - 1, -1))

        kwargs = {"q": q, "axis": axis, "method": method}

        result = apply_ufunc(
            _wrapper,
            self,
            input_core_dims=[dim],
            exclude_dims=set(dim),
            output_core_dims=[["quantile"]],
            output_dtypes=[np.float64],
            dask_gufunc_kwargs=dict(output_sizes={"quantile": len(q)}),
            dask="allowed" if module_available("dask", "2024.11.0") else "parallelized",
            kwargs=kwargs,
        )

        # for backward compatibility
        result = result.transpose("quantile", ...)
        if scalar:
            result = result.squeeze("quantile")
        if keep_attrs:
            result.attrs = self._attrs
        return result

    def rank(self, dim, pct=False):
        """Ranks the data.

        Equal values are assigned a rank that is the average of the ranks that
        would have been otherwise assigned to all of the values within that
        set.  Ranks begin at 1, not 0. If `pct`, computes percentage ranks.

        NaNs in the input array are returned as NaNs.

        The `bottleneck` library is required.

        Parameters
        ----------
        dim : str
            Dimension over which to compute rank.
        pct : bool, optional
            If True, compute percentage ranks, otherwise compute integer ranks.

        Returns
        -------
        ranked : Variable

        See Also
        --------
        Dataset.rank, DataArray.rank
        """
        # This could / should arguably be implemented at the DataArray & Dataset level
        if not OPTIONS["use_bottleneck"]:
            raise RuntimeError(
                "rank requires bottleneck to be enabled."
                " Call `xr.set_options(use_bottleneck=True)` to enable it."
            )

        import bottleneck as bn

        func = bn.nanrankdata if self.dtype.kind == "f" else bn.rankdata
        ranked = xr.apply_ufunc(
            func,
            self,
            input_core_dims=[[dim]],
            output_core_dims=[[dim]],
            dask="parallelized",
            kwargs=dict(axis=-1),
        ).transpose(*self.dims)

        if pct:
            count = self.notnull().sum(dim)
            ranked /= count
        return ranked

    @_deprecate_positional_args("v2024.11.0")
    def rolling_window(
        self,
        dim,
        window,
        window_dim,
        *,
        center=False,
        fill_value=dtypes.NA,
        **kwargs,
    ):
        """
        Make a rolling_window along dim and add a new_dim to the last place.

        Parameters
        ----------
        dim : str
            Dimension over which to compute rolling_window.
            For nd-rolling, should be list of dimensions.
        window : int
            Window size of the rolling
            For nd-rolling, should be list of integers.
        window_dim : str
            New name of the window dimension.
            For nd-rolling, should be list of strings.
        center : bool, default: False
            If True, pad fill_value for both ends. Otherwise, pad in the head
            of the axis.
        fill_value
            value to be filled.
        **kwargs
            Keyword arguments that should be passed to the underlying array type's
            ``sliding_window_view`` function.

        Returns
        -------
        Variable that is a view of the original array with a added dimension of
        size w.
        The return dim: self.dims + (window_dim, )
        The return shape: self.shape + (window, )

        See Also
        --------
        numpy.lib.stride_tricks.sliding_window_view
        dask.array.lib.stride_tricks.sliding_window_view

        Examples
        --------
        >>> v = Variable(("a", "b"), np.arange(8).reshape((2, 4)))
        >>> v.rolling_window("b", 3, "window_dim")
        <xarray.Variable (a: 2, b: 4, window_dim: 3)> Size: 192B
        array([[[nan, nan,  0.],
                [nan,  0.,  1.],
                [ 0.,  1.,  2.],
                [ 1.,  2.,  3.]],
        <BLANKLINE>
               [[nan, nan,  4.],
                [nan,  4.,  5.],
                [ 4.,  5.,  6.],
                [ 5.,  6.,  7.]]])

        >>> v.rolling_window("b", 3, "window_dim", center=True)
        <xarray.Variable (a: 2, b: 4, window_dim: 3)> Size: 192B
        array([[[nan,  0.,  1.],
                [ 0.,  1.,  2.],
                [ 1.,  2.,  3.],
                [ 2.,  3., nan]],
        <BLANKLINE>
               [[nan,  4.,  5.],
                [ 4.,  5.,  6.],
                [ 5.,  6.,  7.],
                [ 6.,  7., nan]]])
        """
        if fill_value is dtypes.NA:  # np.nan is passed
            dtype, fill_value = dtypes.maybe_promote(self.dtype)
            var = duck_array_ops.astype(self, dtype, copy=False)
        else:
            dtype = self.dtype
            var = self

        if utils.is_scalar(dim):
            for name, arg in zip(
                ["window", "window_dim", "center"],
                [window, window_dim, center],
                strict=True,
            ):
                if not utils.is_scalar(arg):
                    raise ValueError(
                        f"Expected {name}={arg!r} to be a scalar like 'dim'."
                    )
            dim = (dim,)

        # dim is now a list
        nroll = len(dim)
        if utils.is_scalar(window):
            window = [window] * nroll
        if utils.is_scalar(window_dim):
            window_dim = [window_dim] * nroll
        if utils.is_scalar(center):
            center = [center] * nroll
        if (
            len(dim) != len(window)
            or len(dim) != len(window_dim)
            or len(dim) != len(center)
        ):
            raise ValueError(
                "'dim', 'window', 'window_dim', and 'center' must be the same length. "
                f"Received dim={dim!r}, window={window!r}, window_dim={window_dim!r},"
                f" and center={center!r}."
            )

        pads = {}
        for d, win, cent in zip(dim, window, center, strict=True):
            if cent:
                start = win // 2  # 10 -> 5,  9 -> 4
                end = win - 1 - start
                pads[d] = (start, end)
            else:
                pads[d] = (win - 1, 0)

        padded = var.pad(pads, mode="constant", constant_values=fill_value)
        axis = self.get_axis_num(dim)
        new_dims = self.dims + tuple(window_dim)
        return Variable(
            new_dims,
            duck_array_ops.sliding_window_view(
                padded.data, window_shape=window, axis=axis, **kwargs
            ),
        )

    def coarsen(
        self, windows, func, boundary="exact", side="left", keep_attrs=None, **kwargs
    ):
        """
        Apply reduction function.
        """
        windows = {k: v for k, v in windows.items() if k in self.dims}

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)

        if keep_attrs:
            _attrs = self.attrs
        else:
            _attrs = None

        if not windows:
            return self._replace(attrs=_attrs)

        reshaped, axes = self.coarsen_reshape(windows, boundary, side)
        if isinstance(func, str):
            name = func
            func = getattr(duck_array_ops, name, None)
            if func is None:
                raise NameError(f"{name} is not a valid method.")

        return self._replace(data=func(reshaped, axis=axes, **kwargs), attrs=_attrs)

    def coarsen_reshape(self, windows, boundary, side):
        """
        Construct a reshaped-array for coarsen
        """
        if not is_dict_like(boundary):
            boundary = dict.fromkeys(windows.keys(), boundary)

        if not is_dict_like(side):
            side = dict.fromkeys(windows.keys(), side)

        # remove unrelated dimensions
        boundary = {k: v for k, v in boundary.items() if k in windows}
        side = {k: v for k, v in side.items() if k in windows}

        for d, window in windows.items():
            if window <= 0:
                raise ValueError(
                    f"window must be > 0. Given {window} for dimension {d}"
                )

        variable = self
        for d, window in windows.items():
            # trim or pad the object
            size = variable.shape[self._get_axis_num(d)]
            n = int(size / window)
            if boundary[d] == "exact":
                if n * window != size:
                    raise ValueError(
                        f"Could not coarsen a dimension of size {size} with "
                        f"window {window} and boundary='exact'. Try a different 'boundary' option."
                    )
            elif boundary[d] == "trim":
                if side[d] == "left":
                    variable = variable.isel({d: slice(0, window * n)})
                else:
                    excess = size - window * n
                    variable = variable.isel({d: slice(excess, None)})
            elif boundary[d] == "pad":  # pad
                pad = window * n - size
                if pad < 0:
                    pad += window
                if side[d] == "left":
                    pad_width = {d: (0, pad)}
                else:
                    pad_width = {d: (pad, 0)}
                variable = variable.pad(pad_width, mode="constant")
            else:
                raise TypeError(
                    f"{boundary[d]} is invalid for boundary. Valid option is 'exact', "
                    "'trim' and 'pad'"
                )

        shape = []
        axes = []
        axis_count = 0
        for i, d in enumerate(variable.dims):
            if d in windows:
                size = variable.shape[i]
                shape.extend((int(size / windows[d]), windows[d]))
                axis_count += 1
                axes.append(i + axis_count)
            else:
                shape.append(variable.shape[i])

        return duck_array_ops.reshape(variable.data, shape), tuple(axes)

    def isnull(self, keep_attrs: bool | None = None):
        """Test each value in the array for whether it is a missing value.

        Returns
        -------
        isnull : Variable
            Same type and shape as object, but the dtype of the data is bool.

        See Also
        --------
        pandas.isnull

        Examples
        --------
        >>> var = xr.Variable("x", [1, np.nan, 3])
        >>> var
        <xarray.Variable (x: 3)> Size: 24B
        array([ 1., nan,  3.])
        >>> var.isnull()
        <xarray.Variable (x: 3)> Size: 3B
        array([False,  True, False])
        """
        from xarray.computation.apply_ufunc import apply_ufunc

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)

        return apply_ufunc(
            duck_array_ops.isnull,
            self,
            dask="allowed",
            keep_attrs=keep_attrs,
        )

    def notnull(self, keep_attrs: bool | None = None):
        """Test each value in the array for whether it is not a missing value.

        Returns
        -------
        notnull : Variable
            Same type and shape as object, but the dtype of the data is bool.

        See Also
        --------
        pandas.notnull

        Examples
        --------
        >>> var = xr.Variable("x", [1, np.nan, 3])
        >>> var
        <xarray.Variable (x: 3)> Size: 24B
        array([ 1., nan,  3.])
        >>> var.notnull()
        <xarray.Variable (x: 3)> Size: 3B
        array([ True, False,  True])
        """
        from xarray.computation.apply_ufunc import apply_ufunc

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)

        return apply_ufunc(
            duck_array_ops.notnull,
            self,
            dask="allowed",
            keep_attrs=keep_attrs,
        )

    @property
    def imag(self) -> Variable:
        """
        The imaginary part of the variable.

        See Also
        --------
        numpy.ndarray.imag
        """
        return self._new(data=self.data.imag)

    @property
    def real(self) -> Variable:
        """
        The real part of the variable.

        See Also
        --------
        numpy.ndarray.real
        """
        return self._new(data=self.data.real)

    def __array_wrap__(self, obj, context=None, return_scalar=False):
        return Variable(self.dims, obj)

    def _unary_op(self, f, *args, **kwargs):
        keep_attrs = kwargs.pop("keep_attrs", None)
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)
        with np.errstate(all="ignore"):
            result = self.__array_wrap__(f(self.data, *args, **kwargs))
            if keep_attrs:
                result.attrs = self.attrs
            return result

    def _binary_op(self, other, f, reflexive=False):
        if isinstance(other, xr.DataTree | xr.DataArray | xr.Dataset):
            return NotImplemented
        if reflexive and issubclass(type(self), type(other)):
            other_data, self_data, dims = _broadcast_compat_data(other, self)
        else:
            self_data, other_data, dims = _broadcast_compat_data(self, other)
        keep_attrs = _get_keep_attrs(default=False)
        attrs = self._attrs if keep_attrs else None
        with np.errstate(all="ignore"):
            new_data = (
                f(self_data, other_data) if not reflexive else f(other_data, self_data)
            )
        result = Variable(dims, new_data, attrs=attrs)
        return result

    def _inplace_binary_op(self, other, f):
        if isinstance(other, xr.Dataset):
            raise TypeError("cannot add a Dataset to a Variable in-place")
        self_data, other_data, dims = _broadcast_compat_data(self, other)
        if dims != self.dims:
            raise ValueError("dimensions cannot change for in-place operations")
        with np.errstate(all="ignore"):
            self.values = f(self_data, other_data)
        return self

    def _to_numeric(self, offset=None, datetime_unit=None, dtype=float):
        """A (private) method to convert datetime array to numeric dtype
        See duck_array_ops.datetime_to_numeric
        """
        numeric_array = duck_array_ops.datetime_to_numeric(
            self.data, offset, datetime_unit, dtype
        )
        return type(self)(self.dims, numeric_array, self._attrs)

    def _unravel_argminmax(
        self,
        argminmax: str,
        dim: Dims,
        axis: int | None,
        keep_attrs: bool | None,
        skipna: bool | None,
    ) -> Variable | dict[Hashable, Variable]:
        """Apply argmin or argmax over one or more dimensions, returning the result as a
        dict of DataArray that can be passed directly to isel.
        """
        if dim is None and axis is None:
            warnings.warn(
                "Behaviour of argmin/argmax with neither dim nor axis argument will "
                "change to return a dict of indices of each dimension. To get a "
                "single, flat index, please use np.argmin(da.data) or "
                "np.argmax(da.data) instead of da.argmin() or da.argmax().",
                DeprecationWarning,
                stacklevel=3,
            )

        argminmax_func = getattr(duck_array_ops, argminmax)

        if dim is ...:
            # In future, should do this also when (dim is None and axis is None)
            dim = self.dims
        if (
            dim is None
            or axis is not None
            or not isinstance(dim, Sequence)
            or isinstance(dim, str)
        ):
            # Return int index if single dimension is passed, and is not part of a
            # sequence
            return self.reduce(
                argminmax_func, dim=dim, axis=axis, keep_attrs=keep_attrs, skipna=skipna
            )

        # Get a name for the new dimension that does not conflict with any existing
        # dimension
        newdimname = "_unravel_argminmax_dim_0"
        count = 1
        while newdimname in self.dims:
            newdimname = f"_unravel_argminmax_dim_{count}"
            count += 1

        stacked = self.stack({newdimname: dim})

        result_dims = stacked.dims[:-1]
        reduce_shape = tuple(self.sizes[d] for d in dim)

        result_flat_indices = stacked.reduce(argminmax_func, axis=-1, skipna=skipna)

        result_unravelled_indices = duck_array_ops.unravel_index(
            result_flat_indices.data, reduce_shape
        )

        result = {
            d: Variable(dims=result_dims, data=i)
            for d, i in zip(dim, result_unravelled_indices, strict=True)
        }

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=False)
        if keep_attrs:
            for v in result.values():
                v.attrs = self.attrs

        return result

    def argmin(
        self,
        dim: Dims = None,
        axis: int | None = None,
        keep_attrs: bool | None = None,
        skipna: bool | None = None,
    ) -> Variable | dict[Hashable, Variable]:
        """Index or indices of the minimum of the Variable over one or more dimensions.
        If a sequence is passed to 'dim', then result returned as dict of Variables,
        which can be passed directly to isel(). If a single str is passed to 'dim' then
        returns a Variable with dtype int.

        If there are multiple minima, the indices of the first one found will be
        returned.

        Parameters
        ----------
        dim : "...", str, Iterable of Hashable or None, optional
            The dimensions over which to find the minimum. By default, finds minimum over
            all dimensions - for now returning an int for backward compatibility, but
            this is deprecated, in future will return a dict with indices for all
            dimensions; to return a dict with all dimensions now, pass '...'.
        axis : int, optional
            Axis over which to apply `argmin`. Only one of the 'dim' and 'axis' arguments
            can be supplied.
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
        result : Variable or dict of Variable

        See Also
        --------
        DataArray.argmin, DataArray.idxmin
        """
        return self._unravel_argminmax("argmin", dim, axis, keep_attrs, skipna)

    def argmax(
        self,
        dim: Dims = None,
        axis: int | None = None,
        keep_attrs: bool | None = None,
        skipna: bool | None = None,
    ) -> Variable | dict[Hashable, Variable]:
        """Index or indices of the maximum of the Variable over one or more dimensions.
        If a sequence is passed to 'dim', then result returned as dict of Variables,
        which can be passed directly to isel(). If a single str is passed to 'dim' then
        returns a Variable with dtype int.

        If there are multiple maxima, the indices of the first one found will be
        returned.

        Parameters
        ----------
        dim : "...", str, Iterable of Hashable or None, optional
            The dimensions over which to find the maximum. By default, finds maximum over
            all dimensions - for now returning an int for backward compatibility, but
            this is deprecated, in future will return a dict with indices for all
            dimensions; to return a dict with all dimensions now, pass '...'.
        axis : int, optional
            Axis over which to apply `argmin`. Only one of the 'dim' and 'axis' arguments
            can be supplied.
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
        result : Variable or dict of Variable

        See Also
        --------
        DataArray.argmax, DataArray.idxmax
        """
        return self._unravel_argminmax("argmax", dim, axis, keep_attrs, skipna)

    def _as_sparse(self, sparse_format=_default, fill_value=_default) -> Variable:
        """
        Use sparse-array as backend.
        """
        from xarray.namedarray._typing import _default as _default_named

        if sparse_format is _default:
            sparse_format = _default_named

        if fill_value is _default:
            fill_value = _default_named

        out = super()._as_sparse(sparse_format, fill_value)
        return cast("Variable", out)

    def _to_dense(self) -> Variable:
        """
        Change backend from sparse to np.array.
        """
        out = super()._to_dense()
        return cast("Variable", out)

    def chunk(  # type: ignore[override]
        self,
        chunks: T_Chunks = {},  # noqa: B006  # even though it's technically unsafe, it is being used intentionally here (#4667)
        name: str | None = None,
        lock: bool | None = None,
        inline_array: bool | None = None,
        chunked_array_type: str | ChunkManagerEntrypoint[Any] | None = None,
        from_array_kwargs: Any = None,
        **chunks_kwargs: Any,
    ) -> Self:
        """Coerce this array's data into a dask array with the given chunks.

        If this variable is a non-dask array, it will be converted to dask
        array. If it's a dask array, it will be rechunked to the given chunk
        sizes.

        If neither chunks is not provided for one or more dimensions, chunk
        sizes along that dimension will not be updated; non-dask arrays will be
        converted into dask arrays with a single block.

        Parameters
        ----------
        chunks : int, tuple or dict, optional
            Chunk sizes along each dimension, e.g., ``5``, ``(5, 5)`` or
            ``{'x': 5, 'y': 5}``.
        name : str, optional
            Used to generate the name for this array in the internal dask
            graph. Does not need not be unique.
        lock : bool, default: False
            Passed on to :py:func:`dask.array.from_array`, if the array is not
            already as dask array.
        inline_array : bool, default: False
            Passed on to :py:func:`dask.array.from_array`, if the array is not
            already as dask array.
        chunked_array_type: str, optional
            Which chunked array type to coerce this datasets' arrays to.
            Defaults to 'dask' if installed, else whatever is registered via the `ChunkManagerEntrypoint` system.
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
        chunked : xarray.Variable

        See Also
        --------
        Variable.chunks
        Variable.chunksizes
        xarray.unify_chunks
        dask.array.from_array
        """

        if from_array_kwargs is None:
            from_array_kwargs = {}

        # TODO deprecate passing these dask-specific arguments explicitly. In future just pass everything via from_array_kwargs
        _from_array_kwargs = consolidate_dask_from_array_kwargs(
            from_array_kwargs,
            name=name,
            lock=lock,
            inline_array=inline_array,
        )

        return super().chunk(
            chunks=chunks,
            chunked_array_type=chunked_array_type,
            from_array_kwargs=_from_array_kwargs,
            **chunks_kwargs,
        )


class IndexVariable(Variable):
    """Wrapper for accommodating a pandas.Index in an xarray.Variable.

    IndexVariable preserve loaded values in the form of a pandas.Index instead
    of a NumPy array. Hence, their values are immutable and must always be one-
    dimensional.

    They also have a name property, which is the name of their sole dimension
    unless another name is given.
    """

    __slots__ = ()

    # TODO: PandasIndexingAdapter doesn't match the array api:
    _data: PandasIndexingAdapter  # type: ignore[assignment]

    def __init__(self, dims, data, attrs=None, encoding=None, fastpath=False):
        super().__init__(dims, data, attrs, encoding, fastpath)
        if self.ndim != 1:
            raise ValueError(f"{type(self).__name__} objects must be 1-dimensional")

        # Unlike in Variable, always eagerly load values into memory
        if not isinstance(self._data, PandasIndexingAdapter):
            self._data = PandasIndexingAdapter(self._data)

    def __dask_tokenize__(self) -> object:
        from dask.base import normalize_token

        # Don't waste time converting pd.Index to np.ndarray
        return normalize_token(
            (type(self), self._dims, self._data.array, self._attrs or None)
        )

    def load(self):
        # data is already loaded into memory for IndexVariable
        return self

    async def load_async(self):
        # data is already loaded into memory for IndexVariable
        return self

    # https://github.com/python/mypy/issues/1465
    @Variable.data.setter  # type: ignore[attr-defined]
    def data(self, data):
        raise ValueError(
            f"Cannot assign to the .data attribute of dimension coordinate a.k.a IndexVariable {self.name!r}. "
            f"Please use DataArray.assign_coords, Dataset.assign_coords or Dataset.assign as appropriate."
        )

    @Variable.values.setter  # type: ignore[attr-defined]
    def values(self, values):
        raise ValueError(
            f"Cannot assign to the .values attribute of dimension coordinate a.k.a IndexVariable {self.name!r}. "
            f"Please use DataArray.assign_coords, Dataset.assign_coords or Dataset.assign as appropriate."
        )

    def chunk(
        self,
        chunks={},  # noqa: B006  # even though it's unsafe, it is being used intentionally here (#4667)
        name=None,
        lock=False,
        inline_array=False,
        chunked_array_type=None,
        from_array_kwargs=None,
    ):
        # Dummy - do not chunk. This method is invoked e.g. by Dataset.chunk()
        return self.copy(deep=False)

    def _as_sparse(self, sparse_format=_default, fill_value=_default):
        # Dummy
        return self.copy(deep=False)

    def _to_dense(self):
        # Dummy
        return self.copy(deep=False)

    def _finalize_indexing_result(self, dims, data):
        if getattr(data, "ndim", 0) != 1:
            # returns Variable rather than IndexVariable if multi-dimensional
            return Variable(dims, data, self._attrs, self._encoding)
        else:
            return self._replace(dims=dims, data=data)

    def __setitem__(self, key, value):
        raise TypeError(f"{type(self).__name__} values cannot be modified")

    @classmethod
    def concat(
        cls,
        variables,
        dim="concat_dim",
        positions=None,
        shortcut=False,
        combine_attrs="override",
    ):
        """Specialized version of Variable.concat for IndexVariable objects.

        This exists because we want to avoid converting Index objects to NumPy
        arrays, if possible.
        """
        from xarray.structure.merge import merge_attrs

        if not isinstance(dim, str):
            (dim,) = dim.dims

        variables = list(variables)
        first_var = variables[0]

        if any(not isinstance(v, cls) for v in variables):
            raise TypeError(
                "IndexVariable.concat requires that all input "
                "variables be IndexVariable objects"
            )

        indexes = [v._data.array for v in variables]

        if not indexes:
            data = []
        else:
            data = indexes[0].append(indexes[1:])

            if positions is not None:
                indices = nputils.inverse_permutation(np.concatenate(positions))
                data = data.take(indices)

        # keep as str if possible as pandas.Index uses object (converts to numpy array)
        data = maybe_coerce_to_str(data, variables)

        attrs = merge_attrs(
            [var.attrs for var in variables], combine_attrs=combine_attrs
        )
        if not shortcut:
            for var in variables:
                if var.dims != first_var.dims:
                    raise ValueError("inconsistent dimensions")

        return cls(first_var.dims, data, attrs)

    def copy(self, deep: bool = True, data: T_DuckArray | ArrayLike | None = None):
        """Returns a copy of this object.

        `deep` is ignored since data is stored in the form of
        pandas.Index, which is already immutable. Dimensions, attributes
        and encodings are always copied.

        Use `data` to create a new object with the same structure as
        original but entirely new data.

        Parameters
        ----------
        deep : bool, default: True
            Deep is ignored when data is given. Whether the data array is
            loaded into memory and copied onto the new object. Default is True.
        data : array_like, optional
            Data to use in the new object. Must have same shape as original.

        Returns
        -------
        object : Variable
            New object with dimensions, attributes, encodings, and optionally
            data copied from original.
        """
        if data is None:
            ndata = self._data

            if deep:
                ndata = copy.deepcopy(ndata, None)

        else:
            ndata = as_compatible_data(data)
            if self.shape != ndata.shape:  # type: ignore[attr-defined]
                raise ValueError(
                    f"Data shape {ndata.shape} must match shape of object {self.shape}"  # type: ignore[attr-defined]
                )

        attrs = copy.deepcopy(self._attrs) if deep else copy.copy(self._attrs)
        encoding = copy.deepcopy(self._encoding) if deep else copy.copy(self._encoding)

        return self._replace(data=ndata, attrs=attrs, encoding=encoding)

    def equals(self, other, equiv=None):
        # if equiv is specified, super up
        if equiv is not None:
            return super().equals(other, equiv)

        # otherwise use the native index equals, rather than looking at _data
        other = getattr(other, "variable", other)
        try:
            return self.dims == other.dims and self._data_equals(other)
        except (TypeError, AttributeError):
            return False

    def _data_equals(self, other):
        return self._to_index().equals(other._to_index())

    def to_index_variable(self) -> IndexVariable:
        """Return this variable as an xarray.IndexVariable"""
        return self.copy(deep=False)

    to_coord = utils.alias(to_index_variable, "to_coord")

    def _to_index(self) -> pd.Index:
        # n.b. creating a new pandas.Index from an old pandas.Index is
        # basically free as pandas.Index objects are immutable.
        # n.b.2. this method returns the multi-index instance for
        # a pandas multi-index level variable.
        assert self.ndim == 1
        index = self._data.array
        if isinstance(index, pd.MultiIndex):
            # set default names for multi-index unnamed levels so that
            # we can safely rename dimension / coordinate later
            valid_level_names = [
                name or f"{self.dims[0]}_level_{i}"
                for i, name in enumerate(index.names)
            ]
            index = index.set_names(valid_level_names)
        else:
            index = index.set_names(self.name)
        return index

    def to_index(self) -> pd.Index:
        """Convert this variable to a pandas.Index"""
        index = self._to_index()
        level = getattr(self._data, "level", None)
        if level is not None:
            # return multi-index level converted to a single index
            return index.get_level_values(level)
        else:
            return index

    @property
    def level_names(self) -> list[Hashable | None] | None:
        """Return MultiIndex level names or None if this IndexVariable has no
        MultiIndex.
        """
        index = self.to_index()
        if isinstance(index, pd.MultiIndex):
            return list(index.names)
        else:
            return None

    def get_level_variable(self, level):
        """Return a new IndexVariable from a given MultiIndex level."""
        if self.level_names is None:
            raise ValueError(f"IndexVariable {self.name!r} has no MultiIndex")
        index = self.to_index()
        return type(self)(self.dims, index.get_level_values(level))

    @property
    def name(self) -> Hashable:
        return self.dims[0]

    @name.setter
    def name(self, value) -> NoReturn:
        raise AttributeError("cannot modify name of IndexVariable in-place")

    def _inplace_binary_op(self, other, f):
        raise TypeError(
            "Values of an IndexVariable are immutable and can not be modified inplace"
        )


def _unified_dims(variables):
    # validate dimensions
    all_dims = {}
    for var in variables:
        var_dims = var.dims
        _raise_if_any_duplicate_dimensions(var_dims, err_context="Broadcasting")

        for d, s in zip(var_dims, var.shape, strict=True):
            if d not in all_dims:
                all_dims[d] = s
            elif all_dims[d] != s:
                raise ValueError(
                    "operands cannot be broadcast together "
                    f"with mismatched lengths for dimension {d!r}: {(all_dims[d], s)}"
                )
    return all_dims


def _broadcast_compat_variables(*variables):
    """Create broadcast compatible variables, with the same dimensions.

    Unlike the result of broadcast_variables(), some variables may have
    dimensions of size 1 instead of the size of the broadcast dimension.
    """
    dims = tuple(_unified_dims(variables))
    return tuple(var.set_dims(dims) if var.dims != dims else var for var in variables)


def broadcast_variables(*variables: Variable) -> tuple[Variable, ...]:
    """Given any number of variables, return variables with matching dimensions
    and broadcast data.

    The data on the returned variables will be a view of the data on the
    corresponding original arrays, but dimensions will be reordered and
    inserted so that both broadcast arrays have the same dimensions. The new
    dimensions are sorted in order of appearance in the first variable's
    dimensions followed by the second variable's dimensions.
    """
    dims_map = _unified_dims(variables)
    dims_tuple = tuple(dims_map)
    return tuple(
        var.set_dims(dims_map) if var.dims != dims_tuple else var for var in variables
    )


def _broadcast_compat_data(self, other):
    if not OPTIONS["arithmetic_broadcast"] and (
        (isinstance(other, Variable) and self.dims != other.dims)
        or (is_duck_array(other) and self.ndim != other.ndim)
    ):
        raise ValueError(
            "Broadcasting is necessary but automatic broadcasting is disabled via "
            "global option `'arithmetic_broadcast'`. "
            "Use `xr.set_options(arithmetic_broadcast=True)` to enable automatic broadcasting."
        )

    if all(hasattr(other, attr) for attr in ["dims", "data", "shape", "encoding"]):
        # `other` satisfies the necessary Variable API for broadcast_variables
        new_self, new_other = _broadcast_compat_variables(self, other)
        self_data = new_self.data
        other_data = new_other.data
        dims = new_self.dims
    else:
        # rely on numpy broadcasting rules
        self_data = self.data
        other_data = other
        dims = self.dims
    return self_data, other_data, dims


def concat(
    variables,
    dim="concat_dim",
    positions=None,
    shortcut=False,
    combine_attrs="override",
):
    """Concatenate variables along a new or existing dimension.

    Parameters
    ----------
    variables : iterable of Variable
        Arrays to stack together. Each variable is expected to have
        matching dimensions and shape except for along the stacked
        dimension.
    dim : str or DataArray, optional
        Name of the dimension to stack along. This can either be a new
        dimension name, in which case it is added along axis=0, or an
        existing dimension name, in which case the location of the
        dimension is unchanged. Where to insert the new dimension is
        determined by the first variable.
    positions : None or list of array-like, optional
        List of integer arrays which specifies the integer positions to which
        to assign each dataset along the concatenated dimension. If not
        supplied, objects are concatenated in the provided order.
    shortcut : bool, optional
        This option is used internally to speed-up groupby operations.
        If `shortcut` is True, some checks of internal consistency between
        arrays to concatenate are skipped.
    combine_attrs : {"drop", "identical", "no_conflicts", "drop_conflicts", \
                     "override"}, default: "override"
        String indicating how to combine attrs of the objects being merged:

        - "drop": empty attrs on returned Dataset.
        - "identical": all attrs must be the same on every object.
        - "no_conflicts": attrs from all objects are combined, any that have
          the same name must also have the same value.
        - "drop_conflicts": attrs from all objects are combined, any that have
          the same name but different values are dropped.
        - "override": skip comparing and copy attrs from the first dataset to
          the result.

    Returns
    -------
    stacked : Variable
        Concatenated Variable formed by stacking all the supplied variables
        along the given dimension.
    """
    variables = list(variables)
    if all(isinstance(v, IndexVariable) for v in variables):
        return IndexVariable.concat(variables, dim, positions, shortcut, combine_attrs)
    else:
        return Variable.concat(variables, dim, positions, shortcut, combine_attrs)


def calculate_dimensions(variables: Mapping[Any, Variable]) -> dict[Hashable, int]:
    """Calculate the dimensions corresponding to a set of variables.

    Returns dictionary mapping from dimension names to sizes. Raises ValueError
    if any of the dimension sizes conflict.
    """
    dims: dict[Hashable, int] = {}
    last_used = {}
    scalar_vars = {k for k, v in variables.items() if not v.dims}
    for k, var in variables.items():
        for dim, size in zip(var.dims, var.shape, strict=True):
            if dim in scalar_vars:
                raise ValueError(
                    f"dimension {dim!r} already exists as a scalar variable"
                )
            if dim not in dims:
                dims[dim] = size
                last_used[dim] = k
            elif dims[dim] != size:
                raise ValueError(
                    f"conflicting sizes for dimension {dim!r}: "
                    f"length {size} on {k!r} and length {dims[dim]} on {last_used!r}"
                )
    return dims
