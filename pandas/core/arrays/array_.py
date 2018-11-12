import numpy as np

from pandas._libs import lib, tslibs

from pandas.core.dtypes.common import is_extension_array_dtype
from pandas.core.dtypes.dtypes import registry
from pandas.core.dtypes.generic import ABCIndexClass, ABCSeries


def array(data, dtype=None, copy=False):
    """
    Create an array.

    .. versionadded:: 0.24.0

    Parameters
    ----------
    data : Sequence[object]
        A sequence of scalar instances for `dtype`. The underlying
        array will be extracted from a Series or Index object.

    dtype : Union[str, np.dtype, ExtensionDtype], optional
        The dtype to use for the array. This may be a NumPy
        dtype, or an extension type registered with pandas using
        :meth:`pandas.api.extensions.register_extension_dtype`.

        By default, the dtype will be inferred from the data
        with :meth:`numpy.array`.

    copy : bool, default False
        Whether to copy the data.

    Returns
    -------
    Array : Union[ndarray, ExtensionArray]

    Examples
    --------
    If a dtype is not specified, `data` is passed through to
    :meth:`numpy.array`, and an ``ndarray`` is returned.

    >>> pd.array([1, 2])
    array([1, 2])

    Or the NumPy dtype can be specified

    >>> pd.array([1, 2], dtype=np.int32)
    array([1, 2], dtype=int32)

    You can use the string alias for `dtype`

    >>> pd.array(['a', 'b', 'a'], dtype='category')
    [a, b, a]
    Categories (2, object): [a, b]

    Or specify the actual dtype

    >>> pd.array(['a', 'b', 'a'],
    ...          dtype=pd.CategoricalDtype(['a', 'b', 'c'], ordered=True))
    [a, b, a]
    Categories (3, object): [a < b < c]

    Because omitting the `dtype` passes the data through to NumPy,
    a mixture of valid integers and NA will return a floating-point
    NumPy array.

    >>> pd.array([1, 2, np.nan])
    array([ 1.,  2., nan])

    To use pandas' nullable :class:`pandas.arrays.IntegerArray`, specify
    the dtype:

    >>> pd.array([1, 2, np.nan], dtype='Int64')
    IntegerArray([1, 2, nan], dtype='Int64')
    """
    from pandas.core.arrays import period_array

    if isinstance(data, (ABCSeries, ABCIndexClass)):
        data = data._values

    # this returns None for not-found dtypes.
    dtype = dtype or registry.find(dtype)

    if is_extension_array_dtype(dtype):
        cls = dtype.construct_array_type()
        return cls._from_sequence(data, dtype=dtype, copy=copy)

    if dtype is None:
        inferred_dtype = lib.infer_dtype(data)
        if inferred_dtype == 'period':
            try:
                return period_array(data)
            except tslibs.IncompatibleFrequency:
                pass  # we return an array below.

        # TODO(DatetimeArray): handle this type
        # TODO(BooleanArray): handle this type

    return np.array(data, dtype=dtype, copy=copy)
