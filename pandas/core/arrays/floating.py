from __future__ import annotations

import numpy as np

from pandas._libs import (
    lib,
    missing as libmissing,
)
from pandas._typing import DtypeObj
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_object_dtype,
)
from pandas.core.dtypes.dtypes import register_extension_dtype

from pandas.core.arrays.numeric import (
    NumericArray,
    NumericDtype,
)


class FloatingDtype(NumericDtype):
    """
    An ExtensionDtype to hold a single size of floating dtype.

    These specific implementations are subclasses of the non-public
    FloatingDtype. For example we have Float32Dtype to represent float32.

    The attributes name & type are set when these subclasses are created.
    """

    def __repr__(self) -> str:
        return f"{self.name}Dtype()"

    @property
    def _is_numeric(self) -> bool:
        return True

    @classmethod
    def construct_array_type(cls) -> type[FloatingArray]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return FloatingArray

    def _get_common_dtype(self, dtypes: list[DtypeObj]) -> DtypeObj | None:
        # for now only handle other floating types
        if not all(isinstance(t, FloatingDtype) for t in dtypes):
            return None
        np_dtype = np.find_common_type(
            # error: Item "ExtensionDtype" of "Union[Any, ExtensionDtype]" has no
            # attribute "numpy_dtype"
            [t.numpy_dtype for t in dtypes],  # type: ignore[union-attr]
            [],
        )
        if np.issubdtype(np_dtype, np.floating):
            return FLOAT_STR_TO_DTYPE[str(np_dtype)]
        return None


def coerce_to_array(
    values, dtype=None, mask=None, copy: bool = False
) -> tuple[np.ndarray, np.ndarray]:
    """
    Coerce the input values array to numpy arrays with a mask.

    Parameters
    ----------
    values : 1D list-like
    dtype : float dtype
    mask : bool 1D array, optional
    copy : bool, default False
        if True, copy the input

    Returns
    -------
    tuple of (values, mask)
    """
    # if values is floating numpy array, preserve its dtype
    if dtype is None and hasattr(values, "dtype"):
        if is_float_dtype(values.dtype):
            dtype = values.dtype

    if dtype is not None:
        if isinstance(dtype, str) and dtype.startswith("Float"):
            # Avoid DeprecationWarning from NumPy about np.dtype("Float64")
            # https://github.com/numpy/numpy/pull/7476
            dtype = dtype.lower()

        if not issubclass(type(dtype), FloatingDtype):
            try:
                dtype = FLOAT_STR_TO_DTYPE[str(np.dtype(dtype))]
            except KeyError as err:
                raise ValueError(f"invalid dtype specified {dtype}") from err

    if isinstance(values, FloatingArray):
        values, mask = values._data, values._mask
        if dtype is not None:
            values = values.astype(dtype.numpy_dtype, copy=False)

        if copy:
            values = values.copy()
            mask = mask.copy()
        return values, mask

    values = np.array(values, copy=copy)
    if is_object_dtype(values.dtype):
        inferred_type = lib.infer_dtype(values, skipna=True)
        if inferred_type == "empty":
            pass
        elif inferred_type == "boolean":
            raise TypeError(f"{values.dtype} cannot be converted to a FloatingDtype")

    elif is_bool_dtype(values) and is_float_dtype(dtype):
        values = np.array(values, dtype=float, copy=copy)

    elif not (is_integer_dtype(values) or is_float_dtype(values)):
        raise TypeError(f"{values.dtype} cannot be converted to a FloatingDtype")

    if values.ndim != 1:
        raise TypeError("values must be a 1D list-like")

    if mask is None:
        mask = libmissing.is_numeric_na(values)

    else:
        assert len(mask) == len(values)

    if not mask.ndim == 1:
        raise TypeError("mask must be a 1D list-like")

    # infer dtype if needed
    if dtype is None:
        dtype = np.dtype("float64")
    else:
        dtype = dtype.type

    # if we are float, let's make sure that we can
    # safely cast

    # we copy as need to coerce here
    # TODO should this be a safe cast?
    if mask.any():
        values = values.copy()
        values[mask] = np.nan
    values = values.astype(dtype, copy=False)  # , casting="safe")

    return values, mask


class FloatingArray(NumericArray):
    """
    Array of floating (optional missing) values.

    .. versionadded:: 1.2.0

    .. warning::

       FloatingArray is currently experimental, and its API or internal
       implementation may change without warning. Especially the behaviour
       regarding NaN (distinct from NA missing values) is subject to change.

    We represent a FloatingArray with 2 numpy arrays:

    - data: contains a numpy float array of the appropriate dtype
    - mask: a boolean array holding a mask on the data, True is missing

    To construct an FloatingArray from generic array-like input, use
    :func:`pandas.array` with one of the float dtypes (see examples).

    See :ref:`integer_na` for more.

    Parameters
    ----------
    values : numpy.ndarray
        A 1-d float-dtype array.
    mask : numpy.ndarray
        A 1-d boolean-dtype array indicating missing values.
    copy : bool, default False
        Whether to copy the `values` and `mask`.

    Attributes
    ----------
    None

    Methods
    -------
    None

    Returns
    -------
    FloatingArray

    Examples
    --------
    Create an FloatingArray with :func:`pandas.array`:

    >>> pd.array([0.1, None, 0.3], dtype=pd.Float32Dtype())
    <FloatingArray>
    [0.1, <NA>, 0.3]
    Length: 3, dtype: Float32

    String aliases for the dtypes are also available. They are capitalized.

    >>> pd.array([0.1, None, 0.3], dtype="Float32")
    <FloatingArray>
    [0.1, <NA>, 0.3]
    Length: 3, dtype: Float32
    """

    # The value used to fill '_data' to avoid upcasting
    _internal_fill_value = 0.0
    # Fill values used for any/all
    _truthy_value = 1.0
    _falsey_value = 0.0

    @cache_readonly
    def dtype(self) -> FloatingDtype:
        return FLOAT_STR_TO_DTYPE[str(self._data.dtype)]

    def __init__(self, values: np.ndarray, mask: np.ndarray, copy: bool = False):
        if not (isinstance(values, np.ndarray) and values.dtype.kind == "f"):
            raise TypeError(
                "values should be floating numpy array. Use "
                "the 'pd.array' function instead"
            )
        if values.dtype == np.float16:
            # If we don't raise here, then accessing self.dtype would raise
            raise TypeError("FloatingArray does not support np.float16 dtype.")

        super().__init__(values, mask, copy=copy)

    @classmethod
    def _coerce_to_array(
        cls, value, *, dtype: DtypeObj, copy: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        return coerce_to_array(value, dtype=dtype, copy=copy)


_dtype_docstring = """
An ExtensionDtype for {dtype} data.

This dtype uses ``pd.NA`` as missing value indicator.

Attributes
----------
None

Methods
-------
None
"""

# create the Dtype


@register_extension_dtype
class Float32Dtype(FloatingDtype):
    type = np.float32
    name = "Float32"
    __doc__ = _dtype_docstring.format(dtype="float32")


@register_extension_dtype
class Float64Dtype(FloatingDtype):
    type = np.float64
    name = "Float64"
    __doc__ = _dtype_docstring.format(dtype="float64")


FLOAT_STR_TO_DTYPE = {
    "float32": Float32Dtype(),
    "float64": Float64Dtype(),
}
