from __future__ import annotations

from typing import Callable

import numpy as np

from pandas._typing import Dtype
from pandas.util._decorators import (
    cache_readonly,
    doc,
)

from pandas.core.dtypes.common import (
    is_dtype_equal,
    is_integer_dtype,
    is_numeric_dtype,
    is_scalar,
    pandas_dtype,
)
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.construction import sanitize_array
from pandas.core.indexes.base import (
    Index,
    maybe_extract_name,
)


class NumericIndex(Index):
    """
    Immutable numeric sequence used for indexing and alignment.

    The basic object storing axis labels for all pandas objects.
    NumericIndex is a special case of `Index` with purely numpy int/uint/float labels.

    .. versionadded:: 1.4.0

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: None)
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
    An NumericIndex instance can **only** contain numpy int64/32/16/8, uint64/32/16/8 or
    float64/32 dtype. In particular, ``NumericIndex`` *can not* hold numpy float16
    dtype or Pandas numeric dtypes (:class:`Int64Dtype`, :class:`Int32Dtype` etc.).
    """

    _typ = "numericindex"
    _values: np.ndarray
    _default_dtype: np.dtype | None = None
    _dtype_validation_metadata: tuple[Callable[..., bool], str] = (
        is_numeric_dtype,
        "numeric type",
    )
    _can_hold_strings = False

    def __new__(
        cls, data=None, dtype: Dtype | None = None, copy: bool = False, name=None
    ) -> NumericIndex:
        name = maybe_extract_name(name, data, cls)

        subarr = cls._ensure_array(data, dtype, copy)
        return cls._simple_new(subarr, name=name)

    @classmethod
    def _ensure_array(cls, data, dtype, copy: bool):
        """
        Ensure we have a valid array to pass to _simple_new.
        """
        cls._validate_dtype(dtype)
        if dtype == np.float16:

            # float16 not supported (no indexing engine)
            raise NotImplementedError("float16 indexes are not supported")

        if not isinstance(data, (np.ndarray, Index)):
            # Coerce to ndarray if not already ndarray or Index
            if is_scalar(data):
                cls._raise_scalar_data_error(data)

            # other iterable of some kind
            if not isinstance(data, (ABCSeries, list, tuple)):
                data = list(data)

            if isinstance(data, (list, tuple)):
                if len(data):
                    data = sanitize_array(data, index=None)
                else:
                    data = np.array([], dtype=np.int64)

        dtype = cls._ensure_dtype(dtype)

        if copy or not is_dtype_equal(data.dtype, dtype):
            # TODO: the try/except below is because it's difficult to predict the error
            # and/or error message from different combinations of data and dtype.
            # Efforts to avoid this try/except welcome.
            # See https://github.com/pandas-dev/pandas/pull/41153#discussion_r676206222
            try:
                subarr = np.array(data, dtype=dtype, copy=copy)
                cls._validate_dtype(subarr.dtype)
            except (TypeError, ValueError):
                raise ValueError(f"data is not compatible with {cls.__name__}")
            cls._assert_safe_casting(data, subarr)
        else:
            subarr = data

        if subarr.ndim > 1:
            # GH#13601, GH#20285, GH#27125
            raise ValueError("Index data must be 1-dimensional")

        subarr = np.asarray(subarr)
        if subarr.dtype == "float16":
            # float16 not supported (no indexing engine)
            raise NotImplementedError("float16 indexes are not implemented")

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

    @classmethod
    def _ensure_dtype(cls, dtype: Dtype | None) -> np.dtype | None:
        """
        Assumes dtype has already been validated.
        """
        if dtype is None:
            return cls._default_dtype

        dtype = pandas_dtype(dtype)
        if not isinstance(dtype, np.dtype):
            raise TypeError(f"{dtype} not a numpy type")
        elif dtype == np.float16:
            # float16 not supported (no indexing engine)
            raise NotImplementedError("float16 indexes are not supported")

        return dtype

    # ----------------------------------------------------------------
    # Indexing Methods

    @cache_readonly
    @doc(Index._should_fallback_to_positional)
    def _should_fallback_to_positional(self) -> bool:
        return False

    # ----------------------------------------------------------------

    @classmethod
    def _assert_safe_casting(cls, data: np.ndarray, subarr: np.ndarray) -> None:
        """
        Ensure incoming data can be represented with matching signed-ness.

        Needed if the process of casting data from some accepted dtype to the internal
        dtype(s) bears the risk of truncation (e.g. float to int).
        """
        if is_integer_dtype(subarr.dtype):
            if not np.array_equal(data, subarr):
                raise TypeError("Unsafe NumPy casting, you must explicitly cast")
