from __future__ import annotations

from typing import (
    Callable,
    Hashable,
)
import warnings

import numpy as np

from pandas._libs import (
    index as libindex,
    lib,
)
from pandas._typing import (
    Dtype,
    DtypeObj,
)
from pandas.util._decorators import doc

from pandas.core.dtypes.cast import astype_nansafe
from pandas.core.dtypes.common import (
    is_dtype_equal,
    is_extension_array_dtype,
    is_float,
    is_float_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_scalar,
    is_signed_integer_dtype,
    is_unsigned_integer_dtype,
    needs_i8_conversion,
    pandas_dtype,
)
from pandas.core.dtypes.generic import ABCSeries

import pandas.core.common as com
from pandas.core.indexes.base import (
    Index,
    maybe_extract_name,
)

_num_index_shared_docs = {}


class NumericIndex(Index):
    """
    Provide numeric type operations.

    This is an abstract class.
    """

    _values: np.ndarray
    _default_dtype: np.dtype
    _dtype_validation_metadata: tuple[Callable[..., bool], str]

    _is_numeric_dtype = True
    _can_hold_strings = False

    def __new__(cls, data=None, dtype: Dtype | None = None, copy=False, name=None):
        name = maybe_extract_name(name, data, cls)

        subarr = cls._ensure_array(data, dtype, copy)
        return cls._simple_new(subarr, name=name)

    @classmethod
    def _ensure_array(cls, data, dtype, copy: bool):
        """
        Ensure we have a valid array to pass to _simple_new.
        """
        cls._validate_dtype(dtype)

        if not isinstance(data, (np.ndarray, Index)):
            # Coerce to ndarray if not already ndarray or Index
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
        return subarr

    @classmethod
    def _validate_dtype(cls, dtype: Dtype | None) -> None:
        if dtype is None:
            return

        validation_func, expected = cls._dtype_validation_metadata
        if not validation_func(dtype):
            raise ValueError(
                f"Incorrect `dtype` passed: expected {expected}, received {dtype}"
            )

    # ----------------------------------------------------------------
    # Indexing Methods

    @doc(Index._maybe_cast_slice_bound)
    def _maybe_cast_slice_bound(self, label, side: str, kind):
        assert kind in ["loc", "getitem", None]

        # we will try to coerce to integers
        return self._maybe_cast_indexer(label)

    # ----------------------------------------------------------------

    @doc(Index._shallow_copy)
    def _shallow_copy(self, values, name: Hashable = lib.no_default):
        if not self._can_hold_na and values.dtype.kind == "f":
            name = self._name if name is lib.no_default else name
            # Ensure we are not returning an Int64Index with float data:
            return Float64Index._simple_new(values, name=name)
        return super()._shallow_copy(values=values, name=name)

    def _convert_tolerance(self, tolerance, target):
        tolerance = super()._convert_tolerance(tolerance, target)

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

    def _is_comparable_dtype(self, dtype: DtypeObj) -> bool:
        # If we ever have BoolIndex or ComplexIndex, this may need to be tightened
        return is_numeric_dtype(dtype)

    @classmethod
    def _assert_safe_casting(cls, data, subarr):
        """
        Subclasses need to override this only if the process of casting data
        from some accepted dtype to the internal dtype(s) bears the risk of
        truncation (e.g. float to int).
        """
        pass

    @property
    def _is_all_dates(self) -> bool:
        """
        Checks that all the labels are datetime objects.
        """
        return False


_num_index_shared_docs[
    "class_descr"
] = """
    Immutable sequence used for indexing and alignment. The basic object
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

_int64_descr_args = {
    "klass": "Int64Index",
    "ltype": "integer",
    "dtype": "int64",
    "extra": "",
}


class IntegerIndex(NumericIndex):
    """
    This is an abstract class for Int64Index, UInt64Index.
    """

    _default_dtype: np.dtype
    _can_hold_na = False

    @classmethod
    def _assert_safe_casting(cls, data, subarr):
        """
        Ensure incoming data can be represented with matching signed-ness.
        """
        if data.dtype.kind != cls._default_dtype.kind:
            if not np.array_equal(data, subarr):
                raise TypeError("Unsafe NumPy casting, you must explicitly cast")

    def __contains__(self, key) -> bool:
        """
        Check if key is a float and has a decimal. If it has, return False.
        """
        hash(key)
        try:
            if is_float(key) and int(key) != key:
                # otherwise the `key in self._engine` check casts e.g. 1.1 -> 1
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
        warnings.warn(
            "Index.asi8 is deprecated and will be removed in a future version",
            FutureWarning,
            stacklevel=2,
        )
        return self._values.view(self._default_dtype)


class Int64Index(IntegerIndex):
    __doc__ = _num_index_shared_docs["class_descr"] % _int64_descr_args

    _typ = "int64index"
    _engine_type = libindex.Int64Engine
    _default_dtype = np.dtype(np.int64)
    _dtype_validation_metadata = (is_signed_integer_dtype, "signed integer")


_uint64_descr_args = {
    "klass": "UInt64Index",
    "ltype": "unsigned integer",
    "dtype": "uint64",
    "extra": "",
}


class UInt64Index(IntegerIndex):
    __doc__ = _num_index_shared_docs["class_descr"] % _uint64_descr_args

    _typ = "uint64index"
    _engine_type = libindex.UInt64Engine
    _default_dtype = np.dtype(np.uint64)
    _dtype_validation_metadata = (is_unsigned_integer_dtype, "unsigned integer")

    # ----------------------------------------------------------------
    # Indexing Methods

    @doc(Index._convert_arr_indexer)
    def _convert_arr_indexer(self, keyarr):
        # Cast the indexer to uint64 if possible so that the values returned
        # from indexing are also uint64.
        dtype = None
        if is_integer_dtype(keyarr) or (
            lib.infer_dtype(keyarr, skipna=False) == "integer"
        ):
            dtype = np.dtype(np.uint64)

        return com.asarray_tuplesafe(keyarr, dtype=dtype)


_float64_descr_args = {
    "klass": "Float64Index",
    "dtype": "float64",
    "ltype": "float",
    "extra": "",
}


class Float64Index(NumericIndex):
    __doc__ = _num_index_shared_docs["class_descr"] % _float64_descr_args

    _typ = "float64index"
    _engine_type = libindex.Float64Engine
    _default_dtype = np.dtype(np.float64)
    _dtype_validation_metadata = (is_float_dtype, "float")

    @property
    def inferred_type(self) -> str:
        """
        Always 'floating' for ``Float64Index``
        """
        return "floating"

    @doc(Index.astype)
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
            return Int64Index(arr, name=self.name)
        return super().astype(dtype, copy=copy)

    # ----------------------------------------------------------------
    # Indexing Methods

    @doc(Index._should_fallback_to_positional)
    def _should_fallback_to_positional(self) -> bool:
        return False

    @doc(Index._convert_slice_indexer)
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
