from typing import Any

import numpy as np

from pandas._libs import index as libindex, lib
from pandas._typing import Dtype, Label
from pandas.util._decorators import Appender, cache_readonly

from pandas.core.dtypes.cast import astype_nansafe
from pandas.core.dtypes.common import (
    is_bool,
    is_bool_dtype,
    is_dtype_equal,
    is_extension_array_dtype,
    is_float,
    is_float_dtype,
    is_integer_dtype,
    is_scalar,
    is_signed_integer_dtype,
    is_unsigned_integer_dtype,
    needs_i8_conversion,
    pandas_dtype,
)
from pandas.core.dtypes.generic import (
    ABCFloat64Index,
    ABCInt64Index,
    ABCRangeIndex,
    ABCSeries,
    ABCUInt64Index,
)
from pandas.core.dtypes.missing import isna

from pandas.core import algorithms
import pandas.core.common as com
from pandas.core.indexes.base import Index, maybe_extract_name
from pandas.core.ops import get_op_result_name

_num_index_shared_docs = dict()


class NumericIndex(Index):
    """
    Provide numeric type operations.

    This is an abstract class.
    """

    _is_numeric_dtype = True

    def __new__(cls, data=None, dtype=None, copy=False, name=None):
        cls._validate_dtype(dtype)
        name = maybe_extract_name(name, data, cls)

        # Coerce to ndarray if not already ndarray or Index
        if not isinstance(data, (np.ndarray, Index)):
            if is_scalar(data):
                raise cls._scalar_data_error(data)

            # other iterable of some kind
            if not isinstance(data, (ABCSeries, list, tuple)):
                data = list(data)

            data = np.asarray(data, dtype=dtype)

        if issubclass(data.dtype.type, str):
            cls._string_data_error(data)

        if copy or not is_dtype_equal(data.dtype, cls._default_dtype):
            subarr = np.array(data, dtype=cls._default_dtype, copy=copy)
            cls._assert_safe_casting(data, subarr)
        else:
            subarr = data

        if subarr.ndim > 1:
            # GH#13601, GH#20285, GH#27125
            raise ValueError("Index data must be 1-dimensional")

        subarr = np.asarray(subarr)
        return cls._simple_new(subarr, name=name)

    @classmethod
    def _validate_dtype(cls, dtype: Dtype) -> None:
        if dtype is None:
            return
        validation_metadata = {
            "int64index": (is_signed_integer_dtype, "signed integer"),
            "uint64index": (is_unsigned_integer_dtype, "unsigned integer"),
            "float64index": (is_float_dtype, "float"),
            "rangeindex": (is_signed_integer_dtype, "signed integer"),
        }

        validation_func, expected = validation_metadata[cls._typ]
        if not validation_func(dtype):
            raise ValueError(
                f"Incorrect `dtype` passed: expected {expected}, received {dtype}"
            )

    @Appender(Index._maybe_cast_slice_bound.__doc__)
    def _maybe_cast_slice_bound(self, label, side, kind):
        assert kind in ["loc", "getitem", None]

        # we will try to coerce to integers
        return self._maybe_cast_indexer(label)

    @Appender(Index._shallow_copy.__doc__)
    def _shallow_copy(self, values=None, name: Label = lib.no_default):
        if values is not None and not self._can_hold_na and values.dtype.kind == "f":
            name = self.name if name is lib.no_default else name
            # Ensure we are not returning an Int64Index with float data:
            return Float64Index._simple_new(values, name=name)
        return super()._shallow_copy(values=values, name=name)

    def _convert_for_op(self, value):
        """
        Convert value to be insertable to ndarray.
        """
        if is_bool(value) or is_bool_dtype(value):
            # force conversion to object
            # so we don't lose the bools
            raise TypeError

        return value

    def _convert_tolerance(self, tolerance, target):
        tolerance = np.asarray(tolerance)
        if target.size != tolerance.size and tolerance.size > 1:
            raise ValueError("list-like tolerance size must match target index size")
        if not np.issubdtype(tolerance.dtype, np.number):
            if tolerance.ndim > 0:
                raise ValueError(
                    f"tolerance argument for {type(self).__name__} must contain "
                    "numeric elements if it is list type"
                )
            else:
                raise ValueError(
                    f"tolerance argument for {type(self).__name__} must be numeric "
                    f"if it is a scalar: {repr(tolerance)}"
                )
        return tolerance

    @classmethod
    def _assert_safe_casting(cls, data, subarr):
        """
        Subclasses need to override this only if the process of casting data
        from some accepted dtype to the internal dtype(s) bears the risk of
        truncation (e.g. float to int).
        """
        pass

    def _concat_same_dtype(self, indexes, name):
        result = type(indexes[0])(np.concatenate([x._values for x in indexes]))
        return result.rename(name)

    @property
    def is_all_dates(self) -> bool:
        """
        Checks that all the labels are datetime objects.
        """
        return False

    @Appender(Index.insert.__doc__)
    def insert(self, loc: int, item):
        # treat NA values as nans:
        if is_scalar(item) and isna(item):
            item = self._na_value
        return super().insert(loc, item)

    def _union(self, other, sort):
        # Right now, we treat union(int, float) a bit special.
        # See https://github.com/pandas-dev/pandas/issues/26778 for discussion
        # We may change union(int, float) to go to object.
        # float | [u]int -> float  (the special case)
        # <T>   | <T>    -> T
        # <T>   | <U>    -> object
        needs_cast = (is_integer_dtype(self.dtype) and is_float_dtype(other.dtype)) or (
            is_integer_dtype(other.dtype) and is_float_dtype(self.dtype)
        )
        if needs_cast:
            first = self.astype("float")
            second = other.astype("float")
            return first._union(second, sort)
        else:
            return super()._union(other, sort)


_num_index_shared_docs[
    "class_descr"
] = """
    Immutable ndarray implementing an ordered, sliceable set. The basic object
    storing axis labels for all pandas objects. %(klass)s is a special case
    of `Index` with purely %(ltype)s labels. %(extra)s.

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: %(dtype)s)
    copy : bool
        Make a copy of input ndarray.
    name : object
        Name to be stored in the index.

    Attributes
    ----------
    None

    Methods
    -------
    None

    See Also
    --------
    Index : The base pandas Index type.

    Notes
    -----
    An Index instance can **only** contain hashable objects.
"""

_int64_descr_args = dict(klass="Int64Index", ltype="integer", dtype="int64", extra="")


class IntegerIndex(NumericIndex):
    """
    This is an abstract class for Int64Index, UInt64Index.
    """

    _default_dtype: np.dtype

    def __contains__(self, key) -> bool:
        """
        Check if key is a float and has a decimal. If it has, return False.
        """
        hash(key)
        try:
            if is_float(key) and int(key) != key:
                return False
            return key in self._engine
        except (OverflowError, TypeError, ValueError):
            return False

    @property
    def inferred_type(self) -> str:
        """
        Always 'integer' for ``Int64Index`` and ``UInt64Index``
        """
        return "integer"

    @property
    def asi8(self) -> np.ndarray:
        # do not cache or you'll create a memory leak
        return self._values.view(self._default_dtype)


class Int64Index(IntegerIndex):
    __doc__ = _num_index_shared_docs["class_descr"] % _int64_descr_args

    _typ = "int64index"
    _can_hold_na = False
    _engine_type = libindex.Int64Engine
    _default_dtype = np.dtype(np.int64)

    def _wrap_joined_index(self, joined, other):
        name = get_op_result_name(self, other)
        return Int64Index(joined, name=name)

    @classmethod
    def _assert_safe_casting(cls, data, subarr):
        """
        Ensure incoming data can be represented as ints.
        """
        if not issubclass(data.dtype.type, np.signedinteger):
            if not np.array_equal(data, subarr):
                raise TypeError("Unsafe NumPy casting, you must explicitly cast")

    def _is_compatible_with_other(self, other) -> bool:
        return super()._is_compatible_with_other(other) or all(
            isinstance(obj, (ABCInt64Index, ABCFloat64Index, ABCRangeIndex))
            for obj in [self, other]
        )


Int64Index._add_numeric_methods()
Int64Index._add_logical_methods()

_uint64_descr_args = dict(
    klass="UInt64Index", ltype="unsigned integer", dtype="uint64", extra=""
)


class UInt64Index(IntegerIndex):
    __doc__ = _num_index_shared_docs["class_descr"] % _uint64_descr_args

    _typ = "uint64index"
    _can_hold_na = False
    _engine_type = libindex.UInt64Engine
    _default_dtype = np.dtype(np.uint64)

    @Appender(Index._convert_arr_indexer.__doc__)
    def _convert_arr_indexer(self, keyarr):
        # Cast the indexer to uint64 if possible so that the values returned
        # from indexing are also uint64.
        dtype = None
        if is_integer_dtype(keyarr) or (
            lib.infer_dtype(keyarr, skipna=False) == "integer"
        ):
            dtype = np.uint64

        return com.asarray_tuplesafe(keyarr, dtype=dtype)

    @Appender(Index._convert_index_indexer.__doc__)
    def _convert_index_indexer(self, keyarr):
        # Cast the indexer to uint64 if possible so
        # that the values returned from indexing are
        # also uint64.
        if keyarr.is_integer():
            return keyarr.astype(np.uint64)
        return keyarr

    def _wrap_joined_index(self, joined, other):
        name = get_op_result_name(self, other)
        return UInt64Index(joined, name=name)

    @classmethod
    def _assert_safe_casting(cls, data, subarr):
        """
        Ensure incoming data can be represented as uints.
        """
        if not issubclass(data.dtype.type, np.unsignedinteger):
            if not np.array_equal(data, subarr):
                raise TypeError("Unsafe NumPy casting, you must explicitly cast")

    def _is_compatible_with_other(self, other) -> bool:
        return super()._is_compatible_with_other(other) or all(
            isinstance(obj, (ABCUInt64Index, ABCFloat64Index)) for obj in [self, other]
        )


UInt64Index._add_numeric_methods()
UInt64Index._add_logical_methods()

_float64_descr_args = dict(
    klass="Float64Index", dtype="float64", ltype="float", extra=""
)


class Float64Index(NumericIndex):
    __doc__ = _num_index_shared_docs["class_descr"] % _float64_descr_args

    _typ = "float64index"
    _engine_type = libindex.Float64Engine
    _default_dtype = np.float64

    @property
    def inferred_type(self) -> str:
        """
        Always 'floating' for ``Float64Index``
        """
        return "floating"

    @Appender(Index.astype.__doc__)
    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)
        if needs_i8_conversion(dtype):
            raise TypeError(
                f"Cannot convert Float64Index to dtype {dtype}; integer "
                "values are required for conversion"
            )
        elif is_integer_dtype(dtype) and not is_extension_array_dtype(dtype):
            # TODO(jreback); this can change once we have an EA Index type
            # GH 13149
            arr = astype_nansafe(self._values, dtype=dtype)
            return Int64Index(arr)
        return super().astype(dtype, copy=copy)

    # ----------------------------------------------------------------
    # Indexing Methods

    @Appender(Index._should_fallback_to_positional.__doc__)
    def _should_fallback_to_positional(self):
        return False

    @Appender(Index._convert_slice_indexer.__doc__)
    def _convert_slice_indexer(self, key: slice, kind: str):
        assert kind in ["loc", "getitem"]

        # We always treat __getitem__ slicing as label-based
        # translate to locations
        return self.slice_indexer(key.start, key.stop, key.step, kind=kind)

    # ----------------------------------------------------------------

    def _format_native_types(
        self, na_rep="", float_format=None, decimal=".", quoting=None, **kwargs
    ):
        from pandas.io.formats.format import FloatArrayFormatter

        formatter = FloatArrayFormatter(
            self._values,
            na_rep=na_rep,
            float_format=float_format,
            decimal=decimal,
            quoting=quoting,
            fixed_width=False,
        )
        return formatter.get_result_as_array()

    def equals(self, other) -> bool:
        """
        Determines if two Index objects contain the same elements.
        """
        if self is other:
            return True

        if not isinstance(other, Index):
            return False

        # need to compare nans locations and make sure that they are the same
        # since nans don't compare equal this is a bit tricky
        try:
            if not isinstance(other, Float64Index):
                other = self._constructor(other)
            if not is_dtype_equal(self.dtype, other.dtype) or self.shape != other.shape:
                return False
            left, right = self._values, other._values
            return ((left == right) | (self._isnan & other._isnan)).all()
        except (TypeError, ValueError):
            return False

    def __contains__(self, other: Any) -> bool:
        hash(other)
        if super().__contains__(other):
            return True

        return is_float(other) and np.isnan(other) and self.hasnans

    @Appender(Index.get_loc.__doc__)
    def get_loc(self, key, method=None, tolerance=None):
        if is_bool(key):
            # Catch this to avoid accidentally casting to 1.0
            raise KeyError(key)

        if is_float(key) and np.isnan(key):
            nan_idxs = self._nan_idxs
            if not len(nan_idxs):
                raise KeyError(key)
            elif len(nan_idxs) == 1:
                return nan_idxs[0]
            return nan_idxs

        return super().get_loc(key, method=method, tolerance=tolerance)

    @cache_readonly
    def is_unique(self) -> bool:
        return super().is_unique and self._nan_idxs.size < 2

    @Appender(Index.isin.__doc__)
    def isin(self, values, level=None):
        if level is not None:
            self._validate_index_level(level)
        return algorithms.isin(np.array(self), values)

    def _is_compatible_with_other(self, other) -> bool:
        return super()._is_compatible_with_other(other) or all(
            isinstance(
                obj, (ABCInt64Index, ABCFloat64Index, ABCUInt64Index, ABCRangeIndex),
            )
            for obj in [self, other]
        )


Float64Index._add_numeric_methods()
Float64Index._add_logical_methods_disabled()
