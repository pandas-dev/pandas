from __future__ import annotations

import itertools

import numpy as np

from pandas._libs.sparse import IntIndex
from pandas._typing import Dtype

from pandas.core.dtypes.common import (
    is_integer_dtype,
    is_list_like,
    is_object_dtype,
)

from pandas.core.arrays import SparseArray
from pandas.core.arrays.categorical import factorize_from_iterable
from pandas.core.frame import DataFrame
from pandas.core.indexes.api import Index
from pandas.core.series import Series


def get_dummies(
    data,
    prefix=None,
    prefix_sep="_",
    dummy_na: bool = False,
    columns=None,
    sparse: bool = False,
    drop_first: bool = False,
    dtype: Dtype | None = None,
) -> DataFrame:
    """
    Convert categorical variable into dummy/indicator variables.

    Parameters
    ----------
    data : array-like, Series, or DataFrame
        Data of which to get dummy indicators.
    prefix : str, list of str, or dict of str, default None
        String to append DataFrame column names.
        Pass a list with length equal to the number of columns
        when calling get_dummies on a DataFrame. Alternatively, `prefix`
        can be a dictionary mapping column names to prefixes.
    prefix_sep : str, default '_'
        If appending prefix, separator/delimiter to use. Or pass a
        list or dictionary as with `prefix`.
    dummy_na : bool, default False
        Add a column to indicate NaNs, if False NaNs are ignored.
    columns : list-like, default None
        Column names in the DataFrame to be encoded.
        If `columns` is None then all the columns with
        `object`, `string`, or `category` dtype will be converted.
    sparse : bool, default False
        Whether the dummy-encoded columns should be backed by
        a :class:`SparseArray` (True) or a regular NumPy array (False).
    drop_first : bool, default False
        Whether to get k-1 dummies out of k categorical levels by removing the
        first level.
    dtype : dtype, default np.uint8
        Data type for new columns. Only a single dtype is allowed.

    Returns
    -------
    DataFrame
        Dummy-coded data.

    See Also
    --------
    Series.str.get_dummies : Convert Series to dummy codes.

    Notes
    -----
    Reference :ref:`the user guide <reshaping.dummies>` for more examples.

    Examples
    --------
    >>> s = pd.Series(list('abca'))

    >>> pd.get_dummies(s)
       a  b  c
    0  1  0  0
    1  0  1  0
    2  0  0  1
    3  1  0  0

    >>> s1 = ['a', 'b', np.nan]

    >>> pd.get_dummies(s1)
       a  b
    0  1  0
    1  0  1
    2  0  0

    >>> pd.get_dummies(s1, dummy_na=True)
       a  b  NaN
    0  1  0    0
    1  0  1    0
    2  0  0    1

    >>> df = pd.DataFrame({'A': ['a', 'b', 'a'], 'B': ['b', 'a', 'c'],
    ...                    'C': [1, 2, 3]})

    >>> pd.get_dummies(df, prefix=['col1', 'col2'])
       C  col1_a  col1_b  col2_a  col2_b  col2_c
    0  1       1       0       0       1       0
    1  2       0       1       1       0       0
    2  3       1       0       0       0       1

    >>> pd.get_dummies(pd.Series(list('abcaa')))
       a  b  c
    0  1  0  0
    1  0  1  0
    2  0  0  1
    3  1  0  0
    4  1  0  0

    >>> pd.get_dummies(pd.Series(list('abcaa')), drop_first=True)
       b  c
    0  0  0
    1  1  0
    2  0  1
    3  0  0
    4  0  0

    >>> pd.get_dummies(pd.Series(list('abc')), dtype=float)
         a    b    c
    0  1.0  0.0  0.0
    1  0.0  1.0  0.0
    2  0.0  0.0  1.0
    """
    from pandas.core.reshape.concat import concat

    dtypes_to_encode = ["object", "string", "category"]

    if isinstance(data, DataFrame):
        # determine columns being encoded
        if columns is None:
            data_to_encode = data.select_dtypes(include=dtypes_to_encode)
        elif not is_list_like(columns):
            raise TypeError("Input must be a list-like for parameter `columns`")
        else:
            data_to_encode = data[columns]

        # validate prefixes and separator to avoid silently dropping cols
        def check_len(item, name):

            if is_list_like(item):
                if not len(item) == data_to_encode.shape[1]:
                    len_msg = (
                        f"Length of '{name}' ({len(item)}) did not match the "
                        "length of the columns being encoded "
                        f"({data_to_encode.shape[1]})."
                    )
                    raise ValueError(len_msg)

        check_len(prefix, "prefix")
        check_len(prefix_sep, "prefix_sep")

        if isinstance(prefix, str):
            prefix = itertools.cycle([prefix])
        if isinstance(prefix, dict):
            prefix = [prefix[col] for col in data_to_encode.columns]

        if prefix is None:
            prefix = data_to_encode.columns

        # validate separators
        if isinstance(prefix_sep, str):
            prefix_sep = itertools.cycle([prefix_sep])
        elif isinstance(prefix_sep, dict):
            prefix_sep = [prefix_sep[col] for col in data_to_encode.columns]

        with_dummies: list[DataFrame]
        if data_to_encode.shape == data.shape:
            # Encoding the entire df, do not prepend any dropped columns
            with_dummies = []
        elif columns is not None:
            # Encoding only cols specified in columns. Get all cols not in
            # columns to prepend to result.
            with_dummies = [data.drop(columns, axis=1)]
        else:
            # Encoding only object and category dtype columns. Get remaining
            # columns to prepend to result.
            with_dummies = [data.select_dtypes(exclude=dtypes_to_encode)]

        for (col, pre, sep) in zip(data_to_encode.items(), prefix, prefix_sep):
            # col is (column_name, column), use just column data here
            dummy = _get_dummies_1d(
                col[1],
                prefix=pre,
                prefix_sep=sep,
                dummy_na=dummy_na,
                sparse=sparse,
                drop_first=drop_first,
                dtype=dtype,
            )
            with_dummies.append(dummy)
        result = concat(with_dummies, axis=1)
    else:
        result = _get_dummies_1d(
            data,
            prefix,
            prefix_sep,
            dummy_na,
            sparse=sparse,
            drop_first=drop_first,
            dtype=dtype,
        )
    return result


def _get_dummies_1d(
    data,
    prefix,
    prefix_sep="_",
    dummy_na: bool = False,
    sparse: bool = False,
    drop_first: bool = False,
    dtype: Dtype | None = None,
) -> DataFrame:
    from pandas.core.reshape.concat import concat

    # Series avoids inconsistent NaN handling
    codes, levels = factorize_from_iterable(Series(data))

    if dtype is None:
        dtype = np.dtype(np.uint8)
    # error: Argument 1 to "dtype" has incompatible type "Union[ExtensionDtype, str,
    # dtype[Any], Type[object]]"; expected "Type[Any]"
    dtype = np.dtype(dtype)  # type: ignore[arg-type]

    if is_object_dtype(dtype):
        raise ValueError("dtype=object is not a valid dtype for get_dummies")

    def get_empty_frame(data) -> DataFrame:
        index: Index | np.ndarray
        if isinstance(data, Series):
            index = data.index
        else:
            index = Index(range(len(data)))
        return DataFrame(index=index)

    # if all NaN
    if not dummy_na and len(levels) == 0:
        return get_empty_frame(data)

    codes = codes.copy()
    if dummy_na:
        codes[codes == -1] = len(levels)
        levels = levels.insert(len(levels), np.nan)

    # if dummy_na, we just fake a nan level. drop_first will drop it again
    if drop_first and len(levels) == 1:
        return get_empty_frame(data)

    number_of_cols = len(levels)

    if prefix is None:
        dummy_cols = levels
    else:
        dummy_cols = Index([f"{prefix}{prefix_sep}{level}" for level in levels])

    index: Index | None
    if isinstance(data, Series):
        index = data.index
    else:
        index = None

    if sparse:

        fill_value: bool | float | int
        if is_integer_dtype(dtype):
            fill_value = 0
        elif dtype == np.dtype(bool):
            fill_value = False
        else:
            fill_value = 0.0

        sparse_series = []
        N = len(data)
        sp_indices: list[list] = [[] for _ in range(len(dummy_cols))]
        mask = codes != -1
        codes = codes[mask]
        n_idx = np.arange(N)[mask]

        for ndx, code in zip(n_idx, codes):
            sp_indices[code].append(ndx)

        if drop_first:
            # remove first categorical level to avoid perfect collinearity
            # GH12042
            sp_indices = sp_indices[1:]
            dummy_cols = dummy_cols[1:]
        for col, ixs in zip(dummy_cols, sp_indices):
            sarr = SparseArray(
                np.ones(len(ixs), dtype=dtype),
                sparse_index=IntIndex(N, ixs),
                fill_value=fill_value,
                dtype=dtype,
            )
            sparse_series.append(Series(data=sarr, index=index, name=col))

        return concat(sparse_series, axis=1, copy=False)

    else:
        # take on axis=1 + transpose to ensure ndarray layout is column-major
        dummy_mat = np.eye(number_of_cols, dtype=dtype).take(codes, axis=1).T

        if not dummy_na:
            # reset NaN GH4446
            dummy_mat[codes == -1] = 0

        if drop_first:
            # remove first GH12042
            dummy_mat = dummy_mat[:, 1:]
            dummy_cols = dummy_cols[1:]
        return DataFrame(dummy_mat, index=index, columns=dummy_cols)
