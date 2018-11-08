import numpy as np

from pandas.core.dtypes.dtypes import registry
from pandas.core.dtypes.common import is_extension_array_dtype
from pandas.core.dtypes.generic import ABCIndexClass, ABCSeries


def array(data, dtype=None, copy=False):
    """
    Create an array.

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
    >>> pd.array([1, 2])
    array([1, 2])

    >>> pd.array(['a', 'b', 'a'], dtype='category')
    [a, b, a]
    Categories (2, object): [a, b]
    """
    if isinstance(data, (ABCSeries, ABCIndexClass)):
        data = data._values

    # this returns None for not-found dtypes.
    dtype = registry.find(dtype) or dtype

    if is_extension_array_dtype(dtype):
        cls = dtype.construct_array_type()
        return cls._from_sequence(data, dtype=dtype, copy=copy)

    return np.array(data, dtype=dtype, copy=copy)
