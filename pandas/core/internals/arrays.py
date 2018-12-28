"""
Methods for cleaning, validating, and unboxing arrays.
"""
from pandas.core.dtypes.generic import ABCIndexClass, ABCPandasArray, ABCSeries


def extract_array(obj, extract_numpy=False):
    """
    Extract the ndarray or ExtensionArray from a Series or Index.

    For all other types, `obj` is just returned as is.

    Parameters
    ----------
    obj : object
        For Series / Index, the underlying ExtensionArray is unboxed.
        For Numpy-backed ExtensionArrays, the ndarray is extracted.

    extract_numpy : bool, default False
        Whether to extract the ndarray from a PandasArray

    Returns
    -------
    arr : object

    Examples
    --------
    >>> extract_array(pd.Series(['a', 'b', 'c'], dtype='category'))
    [a, b, c]
    Categories (3, object): [a, b, c]

    Other objects like lists, arrays, and DataFrames are just passed through.

    >>> extract_array([1, 2, 3])
    [1, 2, 3]

    For an ndarray-backed Series / Index a PandasArray is returned.

    >>> extract_array(pd.Series([1, 2, 3]))
    <PandasArray>
    [1, 2, 3]
    Length: 3, dtype: int64

    To extract all the way down to the ndarray, pass ``extract_numpy=True``.

    >>> extract_array(pd.Series([1, 2, 3]), extract_numpy=True)
    array([1, 2, 3])
    """
    if isinstance(obj, (ABCIndexClass, ABCSeries)):
        obj = obj.array

    if extract_numpy and isinstance(obj, ABCPandasArray):
        obj = obj.to_numpy()

    return obj
