from __future__ import annotations

from copy import deepcopy
import operator
import re
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Literal,
    Sequence,
    TypeVar,
    cast,
)

import numpy as np

from pandas._libs import lib
from pandas._typing import (
    ArrayLike,
    AxisInt,
    Dtype,
    FillnaOptions,
    Iterator,
    NpDtype,
    PositionalIndexer,
    Scalar,
    SortKind,
    TakeIndexer,
    TimeAmbiguous,
    TimeNonexistent,
    npt,
)
from pandas.compat import (
    pa_version_under7p0,
    pa_version_under8p0,
    pa_version_under9p0,
    pa_version_under11p0,
)
from pandas.util._decorators import doc
from pandas.util._validators import validate_fillna_kwargs

from pandas.core.dtypes.common import (
    is_array_like,
    is_bool_dtype,
    is_integer,
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
    is_scalar,
)
from pandas.core.dtypes.missing import isna

from pandas.core import roperator
from pandas.core.arraylike import OpsMixin
from pandas.core.arrays.base import (
    ExtensionArray,
    ExtensionArraySupportsAnyAll,
)
import pandas.core.common as com
from pandas.core.indexers import (
    check_array_indexer,
    unpack_tuple_and_ellipses,
    validate_indices,
)
from pandas.core.strings.base import BaseStringArrayMethods

from pandas.tseries.frequencies import to_offset

if not pa_version_under7p0:
    import pyarrow as pa
    import pyarrow.compute as pc

    from pandas.core.arrays.arrow._arrow_utils import fallback_performancewarning
    from pandas.core.arrays.arrow.dtype import ArrowDtype

    ARROW_CMP_FUNCS = {
        "eq": pc.equal,
        "ne": pc.not_equal,
        "lt": pc.less,
        "gt": pc.greater,
        "le": pc.less_equal,
        "ge": pc.greater_equal,
    }

    ARROW_LOGICAL_FUNCS = {
        "and_": pc.and_kleene,
        "rand_": lambda x, y: pc.and_kleene(y, x),
        "or_": pc.or_kleene,
        "ror_": lambda x, y: pc.or_kleene(y, x),
        "xor": pc.xor,
        "rxor": lambda x, y: pc.xor(y, x),
    }

    def cast_for_truediv(
        arrow_array: pa.ChunkedArray, pa_object: pa.Array | pa.Scalar
    ) -> pa.ChunkedArray:
        # Ensure int / int -> float mirroring Python/Numpy behavior
        # as pc.divide_checked(int, int) -> int
        if pa.types.is_integer(arrow_array.type) and pa.types.is_integer(
            pa_object.type
        ):
            return arrow_array.cast(pa.float64())
        return arrow_array

    def floordiv_compat(
        left: pa.ChunkedArray | pa.Array | pa.Scalar,
        right: pa.ChunkedArray | pa.Array | pa.Scalar,
    ) -> pa.ChunkedArray:
        # Ensure int // int -> int mirroring Python/Numpy behavior
        # as pc.floor(pc.divide_checked(int, int)) -> float
        result = pc.floor(pc.divide(left, right))
        if pa.types.is_integer(left.type) and pa.types.is_integer(right.type):
            result = result.cast(left.type)
        return result

    ARROW_ARITHMETIC_FUNCS = {
        "add": pc.add_checked,
        "radd": lambda x, y: pc.add_checked(y, x),
        "sub": pc.subtract_checked,
        "rsub": lambda x, y: pc.subtract_checked(y, x),
        "mul": pc.multiply_checked,
        "rmul": lambda x, y: pc.multiply_checked(y, x),
        "truediv": lambda x, y: pc.divide(cast_for_truediv(x, y), y),
        "rtruediv": lambda x, y: pc.divide(y, cast_for_truediv(x, y)),
        "floordiv": lambda x, y: floordiv_compat(x, y),
        "rfloordiv": lambda x, y: floordiv_compat(y, x),
        "mod": NotImplemented,
        "rmod": NotImplemented,
        "divmod": NotImplemented,
        "rdivmod": NotImplemented,
        "pow": pc.power_checked,
        "rpow": lambda x, y: pc.power_checked(y, x),
    }

if TYPE_CHECKING:
    from pandas._typing import (
        NumpySorter,
        NumpyValueArrayLike,
    )

    from pandas import Series

ArrowExtensionArrayT = TypeVar("ArrowExtensionArrayT", bound="ArrowExtensionArray")


def get_unit_from_pa_dtype(pa_dtype):
    # https://github.com/pandas-dev/pandas/pull/50998#discussion_r1100344804
    if pa_version_under11p0:
        unit = str(pa_dtype).split("[", 1)[-1][:-1]
        if unit not in ["s", "ms", "us", "ns"]:
            raise ValueError(pa_dtype)
        return unit
    return pa_dtype.unit


def to_pyarrow_type(
    dtype: ArrowDtype | pa.DataType | Dtype | None,
) -> pa.DataType | None:
    """
    Convert dtype to a pyarrow type instance.
    """
    if isinstance(dtype, ArrowDtype):
        return dtype.pyarrow_dtype
    elif isinstance(dtype, pa.DataType):
        return dtype
    elif dtype:
        try:
            # Accepts python types too
            # Doesn't handle all numpy types
            return pa.from_numpy_dtype(dtype)
        except pa.ArrowNotImplementedError:
            pass
    return None


class ArrowExtensionArray(
    OpsMixin, ExtensionArraySupportsAnyAll, BaseStringArrayMethods
):
    """
    Pandas ExtensionArray backed by a PyArrow ChunkedArray.

    .. warning::

       ArrowExtensionArray is considered experimental. The implementation and
       parts of the API may change without warning.

    Parameters
    ----------
    values : pyarrow.Array or pyarrow.ChunkedArray

    Attributes
    ----------
    None

    Methods
    -------
    None

    Returns
    -------
    ArrowExtensionArray

    Notes
    -----
    Most methods are implemented using `pyarrow compute functions. <https://arrow.apache.org/docs/python/api/compute.html>`__
    Some methods may either raise an exception or raise a ``PerformanceWarning`` if an
    associated compute function is not available based on the installed version of PyArrow.

    Please install the latest version of PyArrow to enable the best functionality and avoid
    potential bugs in prior versions of PyArrow.

    Examples
    --------
    Create an ArrowExtensionArray with :func:`pandas.array`:

    >>> pd.array([1, 1, None], dtype="int64[pyarrow]")
    <ArrowExtensionArray>
    [1, 1, <NA>]
    Length: 3, dtype: int64[pyarrow]
    """  # noqa: E501 (http link too long)

    _data: pa.ChunkedArray
    _dtype: ArrowDtype

    def __init__(self, values: pa.Array | pa.ChunkedArray) -> None:
        if pa_version_under7p0:
            msg = "pyarrow>=7.0.0 is required for PyArrow backed ArrowExtensionArray."
            raise ImportError(msg)
        if isinstance(values, pa.Array):
            self._data = pa.chunked_array([values])
        elif isinstance(values, pa.ChunkedArray):
            self._data = values
        else:
            raise ValueError(
                f"Unsupported type '{type(values)}' for ArrowExtensionArray"
            )
        self._dtype = ArrowDtype(self._data.type)

    @classmethod
    def _from_sequence(cls, scalars, *, dtype: Dtype | None = None, copy: bool = False):
        """
        Construct a new ExtensionArray from a sequence of scalars.
        """
        pa_dtype = to_pyarrow_type(dtype)
        if isinstance(scalars, cls):
            scalars = scalars._data
        elif not isinstance(scalars, (pa.Array, pa.ChunkedArray)):
            if copy and is_array_like(scalars):
                # pa array should not get updated when numpy array is updated
                scalars = deepcopy(scalars)
            try:
                scalars = pa.array(scalars, type=pa_dtype, from_pandas=True)
            except pa.ArrowInvalid:
                # GH50430: let pyarrow infer type, then cast
                scalars = pa.array(scalars, from_pandas=True)
        if pa_dtype:
            scalars = scalars.cast(pa_dtype)
        return cls(scalars)

    @classmethod
    def _from_sequence_of_strings(
        cls, strings, *, dtype: Dtype | None = None, copy: bool = False
    ):
        """
        Construct a new ExtensionArray from a sequence of strings.
        """
        pa_type = to_pyarrow_type(dtype)
        if (
            pa_type is None
            or pa.types.is_binary(pa_type)
            or pa.types.is_string(pa_type)
        ):
            # pa_type is None: Let pa.array infer
            # pa_type is string/binary: scalars already correct type
            scalars = strings
        elif pa.types.is_timestamp(pa_type):
            from pandas.core.tools.datetimes import to_datetime

            scalars = to_datetime(strings, errors="raise")
        elif pa.types.is_date(pa_type):
            from pandas.core.tools.datetimes import to_datetime

            scalars = to_datetime(strings, errors="raise").date
        elif pa.types.is_duration(pa_type):
            from pandas.core.tools.timedeltas import to_timedelta

            scalars = to_timedelta(strings, errors="raise")
            if pa_type.unit != "ns":
                # GH51175: test_from_sequence_of_strings_pa_array
                # attempt to parse as int64 reflecting pyarrow's
                # duration to string casting behavior
                mask = isna(scalars)
                if not isinstance(strings, (pa.Array, pa.ChunkedArray)):
                    strings = pa.array(strings, type=pa.string(), from_pandas=True)
                strings = pc.if_else(mask, None, strings)
                try:
                    scalars = strings.cast(pa.int64())
                except pa.ArrowInvalid:
                    pass
        elif pa.types.is_time(pa_type):
            from pandas.core.tools.times import to_time

            # "coerce" to allow "null times" (None) to not raise
            scalars = to_time(strings, errors="coerce")
        elif pa.types.is_boolean(pa_type):
            from pandas.core.arrays import BooleanArray

            scalars = BooleanArray._from_sequence_of_strings(strings).to_numpy()
        elif (
            pa.types.is_integer(pa_type)
            or pa.types.is_floating(pa_type)
            or pa.types.is_decimal(pa_type)
        ):
            from pandas.core.tools.numeric import to_numeric

            scalars = to_numeric(strings, errors="raise")
        else:
            raise NotImplementedError(
                f"Converting strings to {pa_type} is not implemented."
            )
        return cls._from_sequence(scalars, dtype=pa_type, copy=copy)

    def __getitem__(self, item: PositionalIndexer):
        """Select a subset of self.

        Parameters
        ----------
        item : int, slice, or ndarray
            * int: The position in 'self' to get.
            * slice: A slice object, where 'start', 'stop', and 'step' are
              integers or None
            * ndarray: A 1-d boolean NumPy ndarray the same length as 'self'

        Returns
        -------
        item : scalar or ExtensionArray

        Notes
        -----
        For scalar ``item``, return a scalar value suitable for the array's
        type. This should be an instance of ``self.dtype.type``.
        For slice ``key``, return an instance of ``ExtensionArray``, even
        if the slice is length 0 or 1.
        For a boolean mask, return an instance of ``ExtensionArray``, filtered
        to the values where ``item`` is True.
        """
        item = check_array_indexer(self, item)

        if isinstance(item, np.ndarray):
            if not len(item):
                # Removable once we migrate StringDtype[pyarrow] to ArrowDtype[string]
                if self._dtype.name == "string" and self._dtype.storage == "pyarrow":
                    pa_dtype = pa.string()
                else:
                    pa_dtype = self._dtype.pyarrow_dtype
                return type(self)(pa.chunked_array([], type=pa_dtype))
            elif is_integer_dtype(item.dtype):
                return self.take(item)
            elif is_bool_dtype(item.dtype):
                return type(self)(self._data.filter(item))
            else:
                raise IndexError(
                    "Only integers, slices and integer or "
                    "boolean arrays are valid indices."
                )
        elif isinstance(item, tuple):
            item = unpack_tuple_and_ellipses(item)

        if item is Ellipsis:
            # TODO: should be handled by pyarrow?
            item = slice(None)

        if is_scalar(item) and not is_integer(item):
            # e.g. "foo" or 2.5
            # exception message copied from numpy
            raise IndexError(
                r"only integers, slices (`:`), ellipsis (`...`), numpy.newaxis "
                r"(`None`) and integer or boolean arrays are valid indices"
            )
        # We are not an array indexer, so maybe e.g. a slice or integer
        # indexer. We dispatch to pyarrow.
        value = self._data[item]
        if isinstance(value, pa.ChunkedArray):
            return type(self)(value)
        else:
            scalar = value.as_py()
            if scalar is None:
                return self._dtype.na_value
            else:
                return scalar

    def __iter__(self) -> Iterator[Any]:
        """
        Iterate over elements of the array.
        """
        na_value = self._dtype.na_value
        for value in self._data:
            val = value.as_py()
            if val is None:
                yield na_value
            else:
                yield val

    def __arrow_array__(self, type=None):
        """Convert myself to a pyarrow ChunkedArray."""
        return self._data

    def __array__(self, dtype: NpDtype | None = None) -> np.ndarray:
        """Correctly construct numpy arrays when passed to `np.asarray()`."""
        return self.to_numpy(dtype=dtype)

    def __invert__(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        return type(self)(pc.invert(self._data))

    def __neg__(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        return type(self)(pc.negate_checked(self._data))

    def __pos__(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        return type(self)(self._data)

    def __abs__(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        return type(self)(pc.abs_checked(self._data))

    # GH 42600: __getstate__/__setstate__ not necessary once
    # https://issues.apache.org/jira/browse/ARROW-10739 is addressed
    def __getstate__(self):
        state = self.__dict__.copy()
        state["_data"] = self._data.combine_chunks()
        return state

    def __setstate__(self, state) -> None:
        state["_data"] = pa.chunked_array(state["_data"])
        self.__dict__.update(state)

    def _cmp_method(self, other, op):
        pc_func = ARROW_CMP_FUNCS[op.__name__]
        if isinstance(other, ArrowExtensionArray):
            result = pc_func(self._data, other._data)
        elif isinstance(other, (np.ndarray, list)):
            result = pc_func(self._data, other)
        elif is_scalar(other):
            try:
                result = pc_func(self._data, pa.scalar(other))
            except (pa.lib.ArrowNotImplementedError, pa.lib.ArrowInvalid):
                mask = isna(self) | isna(other)
                valid = ~mask
                result = np.zeros(len(self), dtype="bool")
                result[valid] = op(np.array(self)[valid], other)
                result = pa.array(result, type=pa.bool_())
                result = pc.if_else(valid, result, None)
        else:
            raise NotImplementedError(
                f"{op.__name__} not implemented for {type(other)}"
            )
        return ArrowExtensionArray(result)

    def _evaluate_op_method(self, other, op, arrow_funcs):
        pa_type = self._data.type
        if (pa.types.is_string(pa_type) or pa.types.is_binary(pa_type)) and op in [
            operator.add,
            roperator.radd,
        ]:
            length = self._data.length()

            seps: list[str] | list[bytes]
            if pa.types.is_string(pa_type):
                seps = [""] * length
            else:
                seps = [b""] * length

            if is_scalar(other):
                other = [other] * length
            elif isinstance(other, type(self)):
                other = other._data
            if op is operator.add:
                result = pc.binary_join_element_wise(self._data, other, seps)
            else:
                result = pc.binary_join_element_wise(other, self._data, seps)
            return type(self)(result)

        pc_func = arrow_funcs[op.__name__]
        if pc_func is NotImplemented:
            raise NotImplementedError(f"{op.__name__} not implemented.")
        if isinstance(other, ArrowExtensionArray):
            result = pc_func(self._data, other._data)
        elif isinstance(other, (np.ndarray, list)):
            result = pc_func(self._data, pa.array(other, from_pandas=True))
        elif is_scalar(other):
            if isna(other) and op.__name__ in ARROW_LOGICAL_FUNCS:
                # pyarrow kleene ops require null to be typed
                pa_scalar = pa.scalar(None, type=self._data.type)
            else:
                pa_scalar = pa.scalar(other)
            result = pc_func(self._data, pa_scalar)
        else:
            raise NotImplementedError(
                f"{op.__name__} not implemented for {type(other)}"
            )
        return type(self)(result)

    def _logical_method(self, other, op):
        return self._evaluate_op_method(other, op, ARROW_LOGICAL_FUNCS)

    def _arith_method(self, other, op):
        return self._evaluate_op_method(other, op, ARROW_ARITHMETIC_FUNCS)

    def equals(self, other) -> bool:
        if not isinstance(other, ArrowExtensionArray):
            return False
        # I'm told that pyarrow makes __eq__ behave like pandas' equals;
        #  TODO: is this documented somewhere?
        return self._data == other._data

    @property
    def dtype(self) -> ArrowDtype:
        """
        An instance of 'ExtensionDtype'.
        """
        return self._dtype

    @property
    def nbytes(self) -> int:
        """
        The number of bytes needed to store this object in memory.
        """
        return self._data.nbytes

    def __len__(self) -> int:
        """
        Length of this array.

        Returns
        -------
        length : int
        """
        return len(self._data)

    def __contains__(self, key) -> bool:
        # https://github.com/pandas-dev/pandas/pull/51307#issuecomment-1426372604
        if isna(key) and key is not self.dtype.na_value:
            if self.dtype.kind == "f" and lib.is_float(key) and isna(key):
                return pc.any(pc.is_nan(self._data)).as_py()

            # e.g. date or timestamp types we do not allow None here to match pd.NA
            return False
            # TODO: maybe complex? object?

        return bool(super().__contains__(key))

    @property
    def _hasna(self) -> bool:
        return self._data.null_count > 0

    def isna(self) -> npt.NDArray[np.bool_]:
        """
        Boolean NumPy array indicating if each value is missing.

        This should return a 1-D array the same length as 'self'.
        """
        return self._data.is_null().to_numpy()

    def any(self, *, skipna: bool = True, **kwargs):
        """
        Return whether any element is truthy.

        Returns False unless there is at least one element that is truthy.
        By default, NAs are skipped. If ``skipna=False`` is specified and
        missing values are present, similar :ref:`Kleene logic <boolean.kleene>`
        is used as for logical operations.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA values. If the entire array is NA and `skipna` is
            True, then the result will be False, as for an empty array.
            If `skipna` is False, the result will still be True if there is
            at least one element that is truthy, otherwise NA will be returned
            if there are NA's present.

        Returns
        -------
        bool or :attr:`pandas.NA`

        See Also
        --------
        ArrowExtensionArray.all : Return whether all elements are truthy.

        Examples
        --------
        The result indicates whether any element is truthy (and by default
        skips NAs):

        >>> pd.array([True, False, True], dtype="boolean[pyarrow]").any()
        True
        >>> pd.array([True, False, pd.NA], dtype="boolean[pyarrow]").any()
        True
        >>> pd.array([False, False, pd.NA], dtype="boolean[pyarrow]").any()
        False
        >>> pd.array([], dtype="boolean[pyarrow]").any()
        False
        >>> pd.array([pd.NA], dtype="boolean[pyarrow]").any()
        False
        >>> pd.array([pd.NA], dtype="float64[pyarrow]").any()
        False

        With ``skipna=False``, the result can be NA if this is logically
        required (whether ``pd.NA`` is True or False influences the result):

        >>> pd.array([True, False, pd.NA], dtype="boolean[pyarrow]").any(skipna=False)
        True
        >>> pd.array([1, 0, pd.NA], dtype="boolean[pyarrow]").any(skipna=False)
        True
        >>> pd.array([False, False, pd.NA], dtype="boolean[pyarrow]").any(skipna=False)
        <NA>
        >>> pd.array([0, 0, pd.NA], dtype="boolean[pyarrow]").any(skipna=False)
        <NA>
        """
        return self._reduce("any", skipna=skipna, **kwargs)

    def all(self, *, skipna: bool = True, **kwargs):
        """
        Return whether all elements are truthy.

        Returns True unless there is at least one element that is falsey.
        By default, NAs are skipped. If ``skipna=False`` is specified and
        missing values are present, similar :ref:`Kleene logic <boolean.kleene>`
        is used as for logical operations.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA values. If the entire array is NA and `skipna` is
            True, then the result will be True, as for an empty array.
            If `skipna` is False, the result will still be False if there is
            at least one element that is falsey, otherwise NA will be returned
            if there are NA's present.

        Returns
        -------
        bool or :attr:`pandas.NA`

        See Also
        --------
        ArrowExtensionArray.any : Return whether any element is truthy.

        Examples
        --------
        The result indicates whether all elements are truthy (and by default
        skips NAs):

        >>> pd.array([True, True, pd.NA], dtype="boolean[pyarrow]").all()
        True
        >>> pd.array([1, 1, pd.NA], dtype="boolean[pyarrow]").all()
        True
        >>> pd.array([True, False, pd.NA], dtype="boolean[pyarrow]").all()
        False
        >>> pd.array([], dtype="boolean[pyarrow]").all()
        True
        >>> pd.array([pd.NA], dtype="boolean[pyarrow]").all()
        True
        >>> pd.array([pd.NA], dtype="float64[pyarrow]").all()
        True

        With ``skipna=False``, the result can be NA if this is logically
        required (whether ``pd.NA`` is True or False influences the result):

        >>> pd.array([True, True, pd.NA], dtype="boolean[pyarrow]").all(skipna=False)
        <NA>
        >>> pd.array([1, 1, pd.NA], dtype="boolean[pyarrow]").all(skipna=False)
        <NA>
        >>> pd.array([True, False, pd.NA], dtype="boolean[pyarrow]").all(skipna=False)
        False
        >>> pd.array([1, 0, pd.NA], dtype="boolean[pyarrow]").all(skipna=False)
        False
        """
        return self._reduce("all", skipna=skipna, **kwargs)

    def argsort(
        self,
        *,
        ascending: bool = True,
        kind: SortKind = "quicksort",
        na_position: str = "last",
        **kwargs,
    ) -> np.ndarray:
        order = "ascending" if ascending else "descending"
        null_placement = {"last": "at_end", "first": "at_start"}.get(na_position, None)
        if null_placement is None:
            raise ValueError(f"invalid na_position: {na_position}")

        result = pc.array_sort_indices(
            self._data, order=order, null_placement=null_placement
        )
        np_result = result.to_numpy()
        return np_result.astype(np.intp, copy=False)

    def _argmin_max(self, skipna: bool, method: str) -> int:
        if self._data.length() in (0, self._data.null_count) or (
            self._hasna and not skipna
        ):
            # For empty or all null, pyarrow returns -1 but pandas expects TypeError
            # For skipna=False and data w/ null, pandas expects NotImplementedError
            # let ExtensionArray.arg{max|min} raise
            return getattr(super(), f"arg{method}")(skipna=skipna)

        data = self._data
        if pa.types.is_duration(data.type):
            data = data.cast(pa.int64())

        value = getattr(pc, method)(data, skip_nulls=skipna)
        return pc.index(data, value).as_py()

    def argmin(self, skipna: bool = True) -> int:
        return self._argmin_max(skipna, "min")

    def argmax(self, skipna: bool = True) -> int:
        return self._argmin_max(skipna, "max")

    def copy(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        """
        Return a shallow copy of the array.

        Underlying ChunkedArray is immutable, so a deep copy is unnecessary.

        Returns
        -------
        type(self)
        """
        return type(self)(self._data)

    def dropna(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        """
        Return ArrowExtensionArray without NA values.

        Returns
        -------
        ArrowExtensionArray
        """
        return type(self)(pc.drop_null(self._data))

    @doc(ExtensionArray.fillna)
    def fillna(
        self: ArrowExtensionArrayT,
        value: object | ArrayLike | None = None,
        method: FillnaOptions | None = None,
        limit: int | None = None,
    ) -> ArrowExtensionArrayT:
        value, method = validate_fillna_kwargs(value, method)

        if limit is not None:
            return super().fillna(value=value, method=method, limit=limit)

        if method is not None:
            fallback_performancewarning()
            return super().fillna(value=value, method=method, limit=limit)

        if is_array_like(value):
            value = cast(ArrayLike, value)
            if len(value) != len(self):
                raise ValueError(
                    f"Length of 'value' does not match. Got ({len(value)}) "
                    f" expected {len(self)}"
                )

        def convert_fill_value(value, pa_type, dtype):
            if value is None:
                return value
            if isinstance(value, (pa.Scalar, pa.Array, pa.ChunkedArray)):
                return value
            if is_array_like(value):
                pa_box = pa.array
            else:
                pa_box = pa.scalar
            try:
                value = pa_box(value, type=pa_type, from_pandas=True)
            except pa.ArrowTypeError as err:
                msg = f"Invalid value '{str(value)}' for dtype {dtype}"
                raise TypeError(msg) from err
            return value

        fill_value = convert_fill_value(value, self._data.type, self.dtype)

        try:
            if method is None:
                return type(self)(pc.fill_null(self._data, fill_value=fill_value))
            elif method == "pad":
                return type(self)(pc.fill_null_forward(self._data))
            elif method == "backfill":
                return type(self)(pc.fill_null_backward(self._data))
        except pa.ArrowNotImplementedError:
            # ArrowNotImplementedError: Function 'coalesce' has no kernel
            #   matching input types (duration[ns], duration[ns])
            # TODO: remove try/except wrapper if/when pyarrow implements
            #   a kernel for duration types.
            pass

        return super().fillna(value=value, method=method, limit=limit)

    def isin(self, values) -> npt.NDArray[np.bool_]:
        # short-circuit to return all False array.
        if not len(values):
            return np.zeros(len(self), dtype=bool)

        result = pc.is_in(self._data, value_set=pa.array(values, from_pandas=True))
        # pyarrow 2.0.0 returned nulls, so we explicitly specify dtype to convert nulls
        # to False
        return np.array(result, dtype=np.bool_)

    def _values_for_factorize(self) -> tuple[np.ndarray, Any]:
        """
        Return an array and missing value suitable for factorization.

        Returns
        -------
        values : ndarray
        na_value : pd.NA

        Notes
        -----
        The values returned by this method are also used in
        :func:`pandas.util.hash_pandas_object`.
        """
        values = self._data.to_numpy()
        return values, self.dtype.na_value

    @doc(ExtensionArray.factorize)
    def factorize(
        self,
        use_na_sentinel: bool = True,
    ) -> tuple[np.ndarray, ExtensionArray]:
        null_encoding = "mask" if use_na_sentinel else "encode"

        pa_type = self._data.type
        if pa.types.is_duration(pa_type):
            # https://github.com/apache/arrow/issues/15226#issuecomment-1376578323
            data = self._data.cast(pa.int64())
        else:
            data = self._data

        encoded = data.dictionary_encode(null_encoding=null_encoding)
        if encoded.length() == 0:
            indices = np.array([], dtype=np.intp)
            uniques = type(self)(pa.chunked_array([], type=encoded.type.value_type))
        else:
            pa_indices = encoded.combine_chunks().indices
            if pa_indices.null_count > 0:
                pa_indices = pc.fill_null(pa_indices, -1)
            indices = pa_indices.to_numpy(zero_copy_only=False, writable=True).astype(
                np.intp, copy=False
            )
            uniques = type(self)(encoded.chunk(0).dictionary)

        if pa.types.is_duration(pa_type):
            uniques = cast(ArrowExtensionArray, uniques.astype(self.dtype))
        return indices, uniques

    def reshape(self, *args, **kwargs):
        raise NotImplementedError(
            f"{type(self)} does not support reshape "
            f"as backed by a 1D pyarrow.ChunkedArray."
        )

    def round(
        self: ArrowExtensionArrayT, decimals: int = 0, *args, **kwargs
    ) -> ArrowExtensionArrayT:
        """
        Round each value in the array a to the given number of decimals.

        Parameters
        ----------
        decimals : int, default 0
            Number of decimal places to round to. If decimals is negative,
            it specifies the number of positions to the left of the decimal point.
        *args, **kwargs
            Additional arguments and keywords have no effect.

        Returns
        -------
        ArrowExtensionArray
            Rounded values of the ArrowExtensionArray.

        See Also
        --------
        DataFrame.round : Round values of a DataFrame.
        Series.round : Round values of a Series.
        """
        return type(self)(pc.round(self._data, ndigits=decimals))

    @doc(ExtensionArray.searchsorted)
    def searchsorted(
        self,
        value: NumpyValueArrayLike | ExtensionArray,
        side: Literal["left", "right"] = "left",
        sorter: NumpySorter = None,
    ) -> npt.NDArray[np.intp] | np.intp:
        if self._hasna:
            raise ValueError(
                "searchsorted requires array to be sorted, which is impossible "
                "with NAs present."
            )
        if isinstance(value, ExtensionArray):
            value = value.astype(object)
        # Base class searchsorted would cast to object, which is *much* slower.
        return self.to_numpy().searchsorted(value, side=side, sorter=sorter)

    def take(
        self,
        indices: TakeIndexer,
        allow_fill: bool = False,
        fill_value: Any = None,
    ) -> ArrowExtensionArray:
        """
        Take elements from an array.

        Parameters
        ----------
        indices : sequence of int or one-dimensional np.ndarray of int
            Indices to be taken.
        allow_fill : bool, default False
            How to handle negative values in `indices`.

            * False: negative values in `indices` indicate positional indices
              from the right (the default). This is similar to
              :func:`numpy.take`.

            * True: negative values in `indices` indicate
              missing values. These values are set to `fill_value`. Any other
              other negative values raise a ``ValueError``.

        fill_value : any, optional
            Fill value to use for NA-indices when `allow_fill` is True.
            This may be ``None``, in which case the default NA value for
            the type, ``self.dtype.na_value``, is used.

            For many ExtensionArrays, there will be two representations of
            `fill_value`: a user-facing "boxed" scalar, and a low-level
            physical NA value. `fill_value` should be the user-facing version,
            and the implementation should handle translating that to the
            physical version for processing the take if necessary.

        Returns
        -------
        ExtensionArray

        Raises
        ------
        IndexError
            When the indices are out of bounds for the array.
        ValueError
            When `indices` contains negative values other than ``-1``
            and `allow_fill` is True.

        See Also
        --------
        numpy.take
        api.extensions.take

        Notes
        -----
        ExtensionArray.take is called by ``Series.__getitem__``, ``.loc``,
        ``iloc``, when `indices` is a sequence of values. Additionally,
        it's called by :meth:`Series.reindex`, or any other method
        that causes realignment, with a `fill_value`.
        """
        # TODO: Remove once we got rid of the (indices < 0) check
        if not is_array_like(indices):
            indices_array = np.asanyarray(indices)
        else:
            # error: Incompatible types in assignment (expression has type
            # "Sequence[int]", variable has type "ndarray")
            indices_array = indices  # type: ignore[assignment]

        if len(self._data) == 0 and (indices_array >= 0).any():
            raise IndexError("cannot do a non-empty take")
        if indices_array.size > 0 and indices_array.max() >= len(self._data):
            raise IndexError("out of bounds value in 'indices'.")

        if allow_fill:
            fill_mask = indices_array < 0
            if fill_mask.any():
                validate_indices(indices_array, len(self._data))
                # TODO(ARROW-9433): Treat negative indices as NULL
                indices_array = pa.array(indices_array, mask=fill_mask)
                result = self._data.take(indices_array)
                if isna(fill_value):
                    return type(self)(result)
                # TODO: ArrowNotImplementedError: Function fill_null has no
                # kernel matching input types (array[string], scalar[string])
                result = type(self)(result)
                result[fill_mask] = fill_value
                return result
                # return type(self)(pc.fill_null(result, pa.scalar(fill_value)))
            else:
                # Nothing to fill
                return type(self)(self._data.take(indices))
        else:  # allow_fill=False
            # TODO(ARROW-9432): Treat negative indices as indices from the right.
            if (indices_array < 0).any():
                # Don't modify in-place
                indices_array = np.copy(indices_array)
                indices_array[indices_array < 0] += len(self._data)
            return type(self)(self._data.take(indices_array))

    @doc(ExtensionArray.to_numpy)
    def to_numpy(
        self,
        dtype: npt.DTypeLike | None = None,
        copy: bool = False,
        na_value: object = lib.no_default,
    ) -> np.ndarray:
        if dtype is None and self._hasna:
            dtype = object
        if na_value is lib.no_default:
            na_value = self.dtype.na_value

        pa_type = self._data.type
        if pa.types.is_temporal(pa_type) and not pa.types.is_date(pa_type):
            # temporal types with units and/or timezones currently
            #  require pandas/python scalars to pass all tests
            # TODO: improve performance (this is slow)
            result = np.array(list(self), dtype=dtype)
        elif is_object_dtype(dtype) and self._hasna:
            result = np.empty(len(self), dtype=object)
            mask = ~self.isna()
            result[mask] = np.asarray(self[mask]._data)
        elif self._hasna:
            data = self.copy()
            data[self.isna()] = na_value
            return np.asarray(data._data, dtype=dtype)
        else:
            result = np.asarray(self._data, dtype=dtype)
            if copy:
                result = result.copy()
        if self._hasna:
            result[self.isna()] = na_value
        return result

    def unique(self: ArrowExtensionArrayT) -> ArrowExtensionArrayT:
        """
        Compute the ArrowExtensionArray of unique values.

        Returns
        -------
        ArrowExtensionArray
        """
        pa_type = self._data.type

        if pa.types.is_duration(pa_type):
            # https://github.com/apache/arrow/issues/15226#issuecomment-1376578323
            data = self._data.cast(pa.int64())
        else:
            data = self._data

        pa_result = pc.unique(data)

        if pa.types.is_duration(pa_type):
            pa_result = pa_result.cast(pa_type)

        return type(self)(pa_result)

    def value_counts(self, dropna: bool = True) -> Series:
        """
        Return a Series containing counts of each unique value.

        Parameters
        ----------
        dropna : bool, default True
            Don't include counts of missing values.

        Returns
        -------
        counts : Series

        See Also
        --------
        Series.value_counts
        """
        pa_type = self._data.type
        if pa.types.is_duration(pa_type):
            # https://github.com/apache/arrow/issues/15226#issuecomment-1376578323
            data = self._data.cast(pa.int64())
        else:
            data = self._data

        from pandas import (
            Index,
            Series,
        )

        vc = data.value_counts()

        values = vc.field(0)
        counts = vc.field(1)
        if dropna and data.null_count > 0:
            mask = values.is_valid()
            values = values.filter(mask)
            counts = counts.filter(mask)

        if pa.types.is_duration(pa_type):
            values = values.cast(pa_type)

        counts = ArrowExtensionArray(counts)

        index = Index(type(self)(values))

        return Series(counts, index=index, name="count", copy=False)

    @classmethod
    def _concat_same_type(
        cls: type[ArrowExtensionArrayT], to_concat
    ) -> ArrowExtensionArrayT:
        """
        Concatenate multiple ArrowExtensionArrays.

        Parameters
        ----------
        to_concat : sequence of ArrowExtensionArrays

        Returns
        -------
        ArrowExtensionArray
        """
        chunks = [array for ea in to_concat for array in ea._data.iterchunks()]
        if to_concat[0].dtype == "string":
            # StringDtype has no attrivute pyarrow_dtype
            pa_dtype = pa.string()
        else:
            pa_dtype = to_concat[0].dtype.pyarrow_dtype
        arr = pa.chunked_array(chunks, type=pa_dtype)
        return cls(arr)

    def _accumulate(
        self, name: str, *, skipna: bool = True, **kwargs
    ) -> ArrowExtensionArray | ExtensionArray:
        """
        Return an ExtensionArray performing an accumulation operation.

        The underlying data type might change.

        Parameters
        ----------
        name : str
            Name of the function, supported values are:
            - cummin
            - cummax
            - cumsum
            - cumprod
        skipna : bool, default True
            If True, skip NA values.
        **kwargs
            Additional keyword arguments passed to the accumulation function.
            Currently, there is no supported kwarg.

        Returns
        -------
        array

        Raises
        ------
        NotImplementedError : subclass does not define accumulations
        """
        pyarrow_name = {
            "cumsum": "cumulative_sum_checked",
        }.get(name, name)
        pyarrow_meth = getattr(pc, pyarrow_name, None)
        if pyarrow_meth is None:
            return super()._accumulate(name, skipna=skipna, **kwargs)

        data_to_accum = self._data

        pa_dtype = data_to_accum.type
        if pa.types.is_duration(pa_dtype):
            data_to_accum = data_to_accum.cast(pa.int64())

        result = pyarrow_meth(data_to_accum, skip_nulls=skipna, **kwargs)

        if pa.types.is_duration(pa_dtype):
            result = result.cast(pa_dtype)

        return type(self)(result)

    def _reduce(self, name: str, *, skipna: bool = True, **kwargs):
        """
        Return a scalar result of performing the reduction operation.

        Parameters
        ----------
        name : str
            Name of the function, supported values are:
            { any, all, min, max, sum, mean, median, prod,
            std, var, sem, kurt, skew }.
        skipna : bool, default True
            If True, skip NaN values.
        **kwargs
            Additional keyword arguments passed to the reduction function.
            Currently, `ddof` is the only supported kwarg.

        Returns
        -------
        scalar

        Raises
        ------
        TypeError : subclass does not define reductions
        """
        pa_type = self._data.type

        data_to_reduce = self._data

        if name in ["any", "all"] and (
            pa.types.is_integer(pa_type)
            or pa.types.is_floating(pa_type)
            or pa.types.is_duration(pa_type)
            or pa.types.is_decimal(pa_type)
        ):
            # pyarrow only supports any/all for boolean dtype, we allow
            #  for other dtypes, matching our non-pyarrow behavior

            if pa.types.is_duration(pa_type):
                data_to_cmp = self._data.cast(pa.int64())
            else:
                data_to_cmp = self._data

            not_eq = pc.not_equal(data_to_cmp, 0)
            data_to_reduce = not_eq

        elif name in ["min", "max", "sum"] and pa.types.is_duration(pa_type):
            data_to_reduce = self._data.cast(pa.int64())

        elif name in ["median", "mean", "std", "sem"] and pa.types.is_temporal(pa_type):
            nbits = pa_type.bit_width
            if nbits == 32:
                data_to_reduce = self._data.cast(pa.int32())
            else:
                data_to_reduce = self._data.cast(pa.int64())

        if name == "sem":

            def pyarrow_meth(data, skip_nulls, **kwargs):
                numerator = pc.stddev(data, skip_nulls=skip_nulls, **kwargs)
                denominator = pc.sqrt_checked(pc.count(self._data))
                return pc.divide_checked(numerator, denominator)

        else:
            pyarrow_name = {
                "median": "approximate_median",
                "prod": "product",
                "std": "stddev",
                "var": "variance",
            }.get(name, name)
            # error: Incompatible types in assignment
            # (expression has type "Optional[Any]", variable has type
            # "Callable[[Any, Any, KwArg(Any)], Any]")
            pyarrow_meth = getattr(pc, pyarrow_name, None)  # type: ignore[assignment]
            if pyarrow_meth is None:
                # Let ExtensionArray._reduce raise the TypeError
                return super()._reduce(name, skipna=skipna, **kwargs)

        # GH51624: pyarrow defaults to min_count=1, pandas behavior is min_count=0
        if name in ["any", "all"] and "min_count" not in kwargs:
            kwargs["min_count"] = 0

        try:
            result = pyarrow_meth(data_to_reduce, skip_nulls=skipna, **kwargs)
        except (AttributeError, NotImplementedError, TypeError) as err:
            msg = (
                f"'{type(self).__name__}' with dtype {self.dtype} "
                f"does not support reduction '{name}' with pyarrow "
                f"version {pa.__version__}. '{name}' may be supported by "
                f"upgrading pyarrow."
            )
            raise TypeError(msg) from err
        if pc.is_null(result).as_py():
            return self.dtype.na_value

        if name in ["min", "max", "sum"] and pa.types.is_duration(pa_type):
            result = result.cast(pa_type)
        if name in ["median", "mean"] and pa.types.is_temporal(pa_type):
            result = result.cast(pa_type)
        if name in ["std", "sem"] and pa.types.is_temporal(pa_type):
            result = result.cast(pa.int64())
            if pa.types.is_duration(pa_type):
                result = result.cast(pa_type)
            elif pa.types.is_time(pa_type):
                unit = get_unit_from_pa_dtype(pa_type)
                result = result.cast(pa.duration(unit))
            elif pa.types.is_date(pa_type):
                # go with closest available unit, i.e. "s"
                result = result.cast(pa.duration("s"))
            else:
                # i.e. timestamp
                result = result.cast(pa.duration(pa_type.unit))

        return result.as_py()

    def __setitem__(self, key, value) -> None:
        """Set one or more values inplace.

        Parameters
        ----------
        key : int, ndarray, or slice
            When called from, e.g. ``Series.__setitem__``, ``key`` will be
            one of

            * scalar int
            * ndarray of integers.
            * boolean ndarray
            * slice object

        value : ExtensionDtype.type, Sequence[ExtensionDtype.type], or object
            value or values to be set of ``key``.

        Returns
        -------
        None
        """
        # GH50085: unwrap 1D indexers
        if isinstance(key, tuple) and len(key) == 1:
            key = key[0]

        key = check_array_indexer(self, key)
        value = self._maybe_convert_setitem_value(value)

        if com.is_null_slice(key):
            # fast path (GH50248)
            data = self._if_else(True, value, self._data)

        elif is_integer(key):
            # fast path
            key = cast(int, key)
            n = len(self)
            if key < 0:
                key += n
            if not 0 <= key < n:
                raise IndexError(
                    f"index {key} is out of bounds for axis 0 with size {n}"
                )
            if is_list_like(value):
                raise ValueError("Length of indexer and values mismatch")
            elif isinstance(value, pa.Scalar):
                value = value.as_py()
            chunks = [
                *self._data[:key].chunks,
                pa.array([value], type=self._data.type, from_pandas=True),
                *self._data[key + 1 :].chunks,
            ]
            data = pa.chunked_array(chunks).combine_chunks()

        elif is_bool_dtype(key):
            key = np.asarray(key, dtype=np.bool_)
            data = self._replace_with_mask(self._data, key, value)

        elif is_scalar(value) or isinstance(value, pa.Scalar):
            mask = np.zeros(len(self), dtype=np.bool_)
            mask[key] = True
            data = self._if_else(mask, value, self._data)

        else:
            indices = np.arange(len(self))[key]
            if len(indices) != len(value):
                raise ValueError("Length of indexer and values mismatch")
            if len(indices) == 0:
                return
            argsort = np.argsort(indices)
            indices = indices[argsort]
            value = value.take(argsort)
            mask = np.zeros(len(self), dtype=np.bool_)
            mask[indices] = True
            data = self._replace_with_mask(self._data, mask, value)

        if isinstance(data, pa.Array):
            data = pa.chunked_array([data])
        self._data = data

    def _rank(
        self,
        *,
        axis: AxisInt = 0,
        method: str = "average",
        na_option: str = "keep",
        ascending: bool = True,
        pct: bool = False,
    ):
        """
        See Series.rank.__doc__.
        """
        if pa_version_under9p0 or axis != 0:
            ranked = super()._rank(
                axis=axis,
                method=method,
                na_option=na_option,
                ascending=ascending,
                pct=pct,
            )
            # keep dtypes consistent with the implementation below
            if method == "average" or pct:
                pa_type = pa.float64()
            else:
                pa_type = pa.uint64()
            result = pa.array(ranked, type=pa_type, from_pandas=True)
            return type(self)(result)

        data = self._data.combine_chunks()
        sort_keys = "ascending" if ascending else "descending"
        null_placement = "at_start" if na_option == "top" else "at_end"
        tiebreaker = "min" if method == "average" else method

        result = pc.rank(
            data,
            sort_keys=sort_keys,
            null_placement=null_placement,
            tiebreaker=tiebreaker,
        )

        if na_option == "keep":
            mask = pc.is_null(self._data)
            null = pa.scalar(None, type=result.type)
            result = pc.if_else(mask, null, result)

        if method == "average":
            result_max = pc.rank(
                data,
                sort_keys=sort_keys,
                null_placement=null_placement,
                tiebreaker="max",
            )
            result_max = result_max.cast(pa.float64())
            result_min = result.cast(pa.float64())
            result = pc.divide(pc.add(result_min, result_max), 2)

        if pct:
            if not pa.types.is_floating(result.type):
                result = result.cast(pa.float64())
            if method == "dense":
                divisor = pc.max(result)
            else:
                divisor = pc.count(result)
            result = pc.divide(result, divisor)

        return type(self)(result)

    def _quantile(
        self: ArrowExtensionArrayT, qs: npt.NDArray[np.float64], interpolation: str
    ) -> ArrowExtensionArrayT:
        """
        Compute the quantiles of self for each quantile in `qs`.

        Parameters
        ----------
        qs : np.ndarray[float64]
        interpolation: str

        Returns
        -------
        same type as self
        """
        pa_dtype = self._data.type

        data = self._data
        if pa.types.is_temporal(pa_dtype):
            # https://github.com/apache/arrow/issues/33769 in these cases
            #  we can cast to ints and back
            nbits = pa_dtype.bit_width
            if nbits == 32:
                data = data.cast(pa.int32())
            else:
                data = data.cast(pa.int64())

        result = pc.quantile(data, q=qs, interpolation=interpolation)

        if pa.types.is_temporal(pa_dtype):
            nbits = pa_dtype.bit_width
            if nbits == 32:
                result = result.cast(pa.int32())
            else:
                result = result.cast(pa.int64())
            result = result.cast(pa_dtype)

        return type(self)(result)

    def _mode(self: ArrowExtensionArrayT, dropna: bool = True) -> ArrowExtensionArrayT:
        """
        Returns the mode(s) of the ExtensionArray.

        Always returns `ExtensionArray` even if only one value.

        Parameters
        ----------
        dropna : bool, default True
            Don't consider counts of NA values.

        Returns
        -------
        same type as self
            Sorted, if possible.
        """
        pa_type = self._data.type
        if pa.types.is_temporal(pa_type):
            nbits = pa_type.bit_width
            if nbits == 32:
                data = self._data.cast(pa.int32())
            elif nbits == 64:
                data = self._data.cast(pa.int64())
            else:
                raise NotImplementedError(pa_type)
        else:
            data = self._data

        if dropna:
            data = data.drop_null()

        res = pc.value_counts(data)
        most_common = res.field("values").filter(
            pc.equal(res.field("counts"), pc.max(res.field("counts")))
        )

        if pa.types.is_temporal(pa_type):
            most_common = most_common.cast(pa_type)

        return type(self)(most_common)

    def _maybe_convert_setitem_value(self, value):
        """Maybe convert value to be pyarrow compatible."""
        if value is None:
            return value
        if isinstance(value, (pa.Scalar, pa.Array, pa.ChunkedArray)):
            return value
        if is_list_like(value):
            pa_box = pa.array
        else:
            pa_box = pa.scalar
        try:
            value = pa_box(value, type=self._data.type, from_pandas=True)
        except pa.ArrowTypeError as err:
            msg = f"Invalid value '{str(value)}' for dtype {self.dtype}"
            raise TypeError(msg) from err
        return value

    @classmethod
    def _if_else(
        cls,
        cond: npt.NDArray[np.bool_] | bool,
        left: ArrayLike | Scalar,
        right: ArrayLike | Scalar,
    ):
        """
        Choose values based on a condition.

        Analogous to pyarrow.compute.if_else, with logic
        to fallback to numpy for unsupported types.

        Parameters
        ----------
        cond : npt.NDArray[np.bool_] or bool
        left : ArrayLike | Scalar
        right : ArrayLike | Scalar

        Returns
        -------
        pa.Array
        """
        try:
            return pc.if_else(cond, left, right)
        except pa.ArrowNotImplementedError:
            pass

        def _to_numpy_and_type(value) -> tuple[np.ndarray, pa.DataType | None]:
            if isinstance(value, (pa.Array, pa.ChunkedArray)):
                pa_type = value.type
            elif isinstance(value, pa.Scalar):
                pa_type = value.type
                value = value.as_py()
            else:
                pa_type = None
            return np.array(value, dtype=object), pa_type

        left, left_type = _to_numpy_and_type(left)
        right, right_type = _to_numpy_and_type(right)
        pa_type = left_type or right_type
        result = np.where(cond, left, right)
        return pa.array(result, type=pa_type, from_pandas=True)

    @classmethod
    def _replace_with_mask(
        cls,
        values: pa.Array | pa.ChunkedArray,
        mask: npt.NDArray[np.bool_] | bool,
        replacements: ArrayLike | Scalar,
    ):
        """
        Replace items selected with a mask.

        Analogous to pyarrow.compute.replace_with_mask, with logic
        to fallback to numpy for unsupported types.

        Parameters
        ----------
        values : pa.Array or pa.ChunkedArray
        mask : npt.NDArray[np.bool_] or bool
        replacements : ArrayLike or Scalar
            Replacement value(s)

        Returns
        -------
        pa.Array or pa.ChunkedArray
        """
        if isinstance(replacements, pa.ChunkedArray):
            # replacements must be array or scalar, not ChunkedArray
            replacements = replacements.combine_chunks()
        if pa_version_under8p0:
            # pc.replace_with_mask seems to be a bit unreliable for versions < 8.0:
            #  version <= 7: segfaults with various types
            #  version <= 6: fails to replace nulls
            if isinstance(replacements, pa.Array):
                indices = np.full(len(values), None)
                indices[mask] = np.arange(len(replacements))
                indices = pa.array(indices, type=pa.int64())
                replacements = replacements.take(indices)
            return cls._if_else(mask, replacements, values)
        if isinstance(values, pa.ChunkedArray) and pa.types.is_boolean(values.type):
            # GH#52059 replace_with_mask segfaults for chunked array
            # https://github.com/apache/arrow/issues/34634
            values = values.combine_chunks()
        try:
            return pc.replace_with_mask(values, mask, replacements)
        except pa.ArrowNotImplementedError:
            pass
        if isinstance(replacements, pa.Array):
            replacements = np.array(replacements, dtype=object)
        elif isinstance(replacements, pa.Scalar):
            replacements = replacements.as_py()
        result = np.array(values, dtype=object)
        result[mask] = replacements
        return pa.array(result, type=values.type, from_pandas=True)

    def _str_count(self, pat: str, flags: int = 0):
        if flags:
            raise NotImplementedError(f"count not implemented with {flags=}")
        return type(self)(pc.count_substring_regex(self._data, pat))

    def _str_pad(
        self,
        width: int,
        side: Literal["left", "right", "both"] = "left",
        fillchar: str = " ",
    ):
        if side == "left":
            pa_pad = pc.utf8_lpad
        elif side == "right":
            pa_pad = pc.utf8_rpad
        elif side == "both":
            pa_pad = pc.utf8_center
        else:
            raise ValueError(
                f"Invalid side: {side}. Side must be one of 'left', 'right', 'both'"
            )
        return type(self)(pa_pad(self._data, width=width, padding=fillchar))

    def _str_contains(
        self, pat, case: bool = True, flags: int = 0, na=None, regex: bool = True
    ):
        if flags:
            raise NotImplementedError(f"contains not implemented with {flags=}")

        if regex:
            pa_contains = pc.match_substring_regex
        else:
            pa_contains = pc.match_substring
        result = pa_contains(self._data, pat, ignore_case=not case)
        if not isna(na):
            result = result.fill_null(na)
        return type(self)(result)

    def _str_startswith(self, pat: str, na=None):
        result = pc.starts_with(self._data, pattern=pat)
        if not isna(na):
            result = result.fill_null(na)
        return type(self)(result)

    def _str_endswith(self, pat: str, na=None):
        result = pc.ends_with(self._data, pattern=pat)
        if not isna(na):
            result = result.fill_null(na)
        return type(self)(result)

    def _str_replace(
        self,
        pat: str | re.Pattern,
        repl: str | Callable,
        n: int = -1,
        case: bool = True,
        flags: int = 0,
        regex: bool = True,
    ):
        if isinstance(pat, re.Pattern) or callable(repl) or not case or flags:
            raise NotImplementedError(
                "replace is not supported with a re.Pattern, callable repl, "
                "case=False, or flags!=0"
            )

        func = pc.replace_substring_regex if regex else pc.replace_substring
        result = func(self._data, pattern=pat, replacement=repl, max_replacements=n)
        return type(self)(result)

    def _str_repeat(self, repeats: int | Sequence[int]):
        if not isinstance(repeats, int):
            raise NotImplementedError(
                f"repeat is not implemented when repeats is {type(repeats).__name__}"
            )
        elif pa_version_under7p0:
            raise NotImplementedError("repeat is not implemented for pyarrow < 7")
        else:
            return type(self)(pc.binary_repeat(self._data, repeats))

    def _str_match(
        self, pat: str, case: bool = True, flags: int = 0, na: Scalar | None = None
    ):
        if not pat.startswith("^"):
            pat = f"^{pat}"
        return self._str_contains(pat, case, flags, na, regex=True)

    def _str_fullmatch(
        self, pat, case: bool = True, flags: int = 0, na: Scalar | None = None
    ):
        if not pat.endswith("$") or pat.endswith("//$"):
            pat = f"{pat}$"
        return self._str_match(pat, case, flags, na)

    def _str_find(self, sub: str, start: int = 0, end: int | None = None):
        if start != 0 and end is not None:
            slices = pc.utf8_slice_codeunits(self._data, start, stop=end)
            result = pc.find_substring(slices, sub)
            not_found = pc.equal(result, -1)
            offset_result = pc.add(result, end - start)
            result = pc.if_else(not_found, result, offset_result)
        elif start == 0 and end is None:
            slices = self._data
            result = pc.find_substring(slices, sub)
        else:
            raise NotImplementedError(
                f"find not implemented with {sub=}, {start=}, {end=}"
            )
        return type(self)(result)

    def _str_get(self, i: int):
        lengths = pc.utf8_length(self._data)
        if i >= 0:
            out_of_bounds = pc.greater_equal(i, lengths)
            start = i
            stop = i + 1
            step = 1
        else:
            out_of_bounds = pc.greater(-i, lengths)
            start = i
            stop = i - 1
            step = -1
        not_out_of_bounds = pc.invert(out_of_bounds.fill_null(True))
        selected = pc.utf8_slice_codeunits(
            self._data, start=start, stop=stop, step=step
        )
        result = pa.array([None] * self._data.length(), type=self._data.type)
        result = pc.if_else(not_out_of_bounds, selected, result)
        return type(self)(result)

    def _str_join(self, sep: str):
        return type(self)(pc.binary_join(self._data, sep))

    def _str_partition(self, sep: str, expand: bool):
        raise NotImplementedError(
            "str.partition not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_rpartition(self, sep: str, expand: bool):
        raise NotImplementedError(
            "str.rpartition not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_slice(
        self, start: int | None = None, stop: int | None = None, step: int | None = None
    ):
        if start is None:
            start = 0
        if step is None:
            step = 1
        return type(self)(
            pc.utf8_slice_codeunits(self._data, start=start, stop=stop, step=step)
        )

    def _str_slice_replace(
        self, start: int | None = None, stop: int | None = None, repl: str | None = None
    ):
        if repl is None:
            repl = ""
        if start is None:
            start = 0
        return type(self)(pc.utf8_replace_slice(self._data, start, stop, repl))

    def _str_isalnum(self):
        return type(self)(pc.utf8_is_alnum(self._data))

    def _str_isalpha(self):
        return type(self)(pc.utf8_is_alpha(self._data))

    def _str_isdecimal(self):
        return type(self)(pc.utf8_is_decimal(self._data))

    def _str_isdigit(self):
        return type(self)(pc.utf8_is_digit(self._data))

    def _str_islower(self):
        return type(self)(pc.utf8_is_lower(self._data))

    def _str_isnumeric(self):
        return type(self)(pc.utf8_is_numeric(self._data))

    def _str_isspace(self):
        return type(self)(pc.utf8_is_space(self._data))

    def _str_istitle(self):
        return type(self)(pc.utf8_is_title(self._data))

    def _str_capitalize(self):
        return type(self)(pc.utf8_capitalize(self._data))

    def _str_title(self):
        return type(self)(pc.utf8_title(self._data))

    def _str_isupper(self):
        return type(self)(pc.utf8_is_upper(self._data))

    def _str_swapcase(self):
        return type(self)(pc.utf8_swapcase(self._data))

    def _str_len(self):
        return type(self)(pc.utf8_length(self._data))

    def _str_lower(self):
        return type(self)(pc.utf8_lower(self._data))

    def _str_upper(self):
        return type(self)(pc.utf8_upper(self._data))

    def _str_strip(self, to_strip=None):
        if to_strip is None:
            result = pc.utf8_trim_whitespace(self._data)
        else:
            result = pc.utf8_trim(self._data, characters=to_strip)
        return type(self)(result)

    def _str_lstrip(self, to_strip=None):
        if to_strip is None:
            result = pc.utf8_ltrim_whitespace(self._data)
        else:
            result = pc.utf8_ltrim(self._data, characters=to_strip)
        return type(self)(result)

    def _str_rstrip(self, to_strip=None):
        if to_strip is None:
            result = pc.utf8_rtrim_whitespace(self._data)
        else:
            result = pc.utf8_rtrim(self._data, characters=to_strip)
        return type(self)(result)

    def _str_removeprefix(self, prefix: str):
        raise NotImplementedError(
            "str.removeprefix not supported with pd.ArrowDtype(pa.string())."
        )
        # TODO: Should work once https://github.com/apache/arrow/issues/14991 is fixed
        # starts_with = pc.starts_with(self._data, pattern=prefix)
        # removed = pc.utf8_slice_codeunits(self._data, len(prefix))
        # result = pc.if_else(starts_with, removed, self._data)
        # return type(self)(result)

    def _str_removesuffix(self, suffix: str):
        ends_with = pc.ends_with(self._data, pattern=suffix)
        removed = pc.utf8_slice_codeunits(self._data, 0, stop=-len(suffix))
        result = pc.if_else(ends_with, removed, self._data)
        return type(self)(result)

    def _str_casefold(self):
        raise NotImplementedError(
            "str.casefold not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_encode(self, encoding, errors: str = "strict"):
        raise NotImplementedError(
            "str.encode not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_extract(self, pat: str, flags: int = 0, expand: bool = True):
        raise NotImplementedError(
            "str.extract not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_findall(self, pat, flags: int = 0):
        raise NotImplementedError(
            "str.findall not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_get_dummies(self, sep: str = "|"):
        raise NotImplementedError(
            "str.get_dummies not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_index(self, sub, start: int = 0, end=None):
        raise NotImplementedError(
            "str.index not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_rindex(self, sub, start: int = 0, end=None):
        raise NotImplementedError(
            "str.rindex not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_normalize(self, form):
        raise NotImplementedError(
            "str.normalize not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_rfind(self, sub, start: int = 0, end=None):
        raise NotImplementedError(
            "str.rfind not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_split(
        self, pat=None, n=-1, expand: bool = False, regex: bool | None = None
    ):
        raise NotImplementedError(
            "str.split not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_rsplit(self, pat=None, n=-1):
        raise NotImplementedError(
            "str.rsplit not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_translate(self, table):
        raise NotImplementedError(
            "str.translate not supported with pd.ArrowDtype(pa.string())."
        )

    def _str_wrap(self, width, **kwargs):
        raise NotImplementedError(
            "str.wrap not supported with pd.ArrowDtype(pa.string())."
        )

    @property
    def _dt_year(self):
        return type(self)(pc.year(self._data))

    @property
    def _dt_day(self):
        return type(self)(pc.day(self._data))

    @property
    def _dt_day_of_week(self):
        return type(self)(pc.day_of_week(self._data))

    _dt_dayofweek = _dt_day_of_week
    _dt_weekday = _dt_day_of_week

    @property
    def _dt_day_of_year(self):
        return type(self)(pc.day_of_year(self._data))

    _dt_dayofyear = _dt_day_of_year

    @property
    def _dt_hour(self):
        return type(self)(pc.hour(self._data))

    def _dt_isocalendar(self):
        return type(self)(pc.iso_calendar(self._data))

    @property
    def _dt_is_leap_year(self):
        return type(self)(pc.is_leap_year(self._data))

    @property
    def _dt_microsecond(self):
        return type(self)(pc.microsecond(self._data))

    @property
    def _dt_minute(self):
        return type(self)(pc.minute(self._data))

    @property
    def _dt_month(self):
        return type(self)(pc.month(self._data))

    @property
    def _dt_nanosecond(self):
        return type(self)(pc.nanosecond(self._data))

    @property
    def _dt_quarter(self):
        return type(self)(pc.quarter(self._data))

    @property
    def _dt_second(self):
        return type(self)(pc.second(self._data))

    @property
    def _dt_date(self):
        return type(self)(self._data.cast(pa.date32()))

    @property
    def _dt_time(self):
        unit = (
            self.dtype.pyarrow_dtype.unit
            if self.dtype.pyarrow_dtype.unit in {"us", "ns"}
            else "ns"
        )
        return type(self)(self._data.cast(pa.time64(unit)))

    @property
    def _dt_tz(self):
        return self.dtype.pyarrow_dtype.tz

    def _dt_strftime(self, format: str):
        return type(self)(pc.strftime(self._data, format=format))

    def _round_temporally(
        self,
        method: Literal["ceil", "floor", "round"],
        freq,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ):
        if ambiguous != "raise":
            raise NotImplementedError("ambiguous is not supported.")
        if nonexistent != "raise":
            raise NotImplementedError("nonexistent is not supported.")
        offset = to_offset(freq)
        if offset is None:
            raise ValueError(f"Must specify a valid frequency: {freq}")
        pa_supported_unit = {
            "A": "year",
            "AS": "year",
            "Q": "quarter",
            "QS": "quarter",
            "M": "month",
            "MS": "month",
            "W": "week",
            "D": "day",
            "H": "hour",
            "T": "minute",
            "S": "second",
            "L": "millisecond",
            "U": "microsecond",
            "N": "nanosecond",
        }
        unit = pa_supported_unit.get(offset._prefix, None)
        if unit is None:
            raise ValueError(f"{freq=} is not supported")
        multiple = offset.n
        rounding_method = getattr(pc, f"{method}_temporal")
        return type(self)(rounding_method(self._data, multiple=multiple, unit=unit))

    def _dt_ceil(
        self,
        freq,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ):
        return self._round_temporally("ceil", freq, ambiguous, nonexistent)

    def _dt_floor(
        self,
        freq,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ):
        return self._round_temporally("floor", freq, ambiguous, nonexistent)

    def _dt_round(
        self,
        freq,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ):
        return self._round_temporally("round", freq, ambiguous, nonexistent)

    def _dt_to_pydatetime(self):
        data = self._data.to_pylist()
        if self._dtype.pyarrow_dtype.unit == "ns":
            data = [None if ts is None else ts.to_pydatetime(warn=False) for ts in data]
        return np.array(data, dtype=object)

    def _dt_tz_localize(
        self,
        tz,
        ambiguous: TimeAmbiguous = "raise",
        nonexistent: TimeNonexistent = "raise",
    ):
        if ambiguous != "raise":
            raise NotImplementedError(f"{ambiguous=} is not supported")
        if nonexistent != "raise":
            raise NotImplementedError(f"{nonexistent=} is not supported")
        if tz is None:
            new_type = pa.timestamp(self.dtype.pyarrow_dtype.unit)
            return type(self)(self._data.cast(new_type))
        pa_tz = str(tz)
        return type(self)(
            self._data.cast(pa.timestamp(self.dtype.pyarrow_dtype.unit, pa_tz))
        )
