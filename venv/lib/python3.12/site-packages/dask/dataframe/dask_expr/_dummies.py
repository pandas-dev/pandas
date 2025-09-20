from __future__ import annotations

import pandas as pd

import dask.dataframe.methods as methods
from dask.dataframe.dask_expr._collection import DataFrame, Series, new_collection
from dask.dataframe.dask_expr._expr import Blockwise
from dask.dataframe.utils import has_known_categories
from dask.utils import get_meta_library


def get_dummies(
    data,
    prefix=None,
    prefix_sep="_",
    dummy_na=False,
    columns=None,
    sparse=False,
    drop_first=False,
    dtype=bool,
    **kwargs,
):
    """
    Convert categorical variable into dummy/indicator variables.

    Data must have category dtype to infer result's ``columns``.

    Parameters
    ----------
    data : Series, or DataFrame
        For Series, the dtype must be categorical.
        For DataFrame, at least one column must be categorical.
    prefix : string, list of strings, or dict of strings, default None
        String to append DataFrame column names.
        Pass a list with length equal to the number of columns
        when calling get_dummies on a DataFrame. Alternatively, `prefix`
        can be a dictionary mapping column names to prefixes.
    prefix_sep : string, default '_'
        If appending prefix, separator/delimiter to use. Or pass a
        list or dictionary as with `prefix.`
    dummy_na : bool, default False
        Add a column to indicate NaNs, if False NaNs are ignored.
    columns : list-like, default None
        Column names in the DataFrame to be encoded.
        If `columns` is None then all the columns with
        `category` dtype will be converted.
    sparse : bool, default False
        Whether the dummy columns should be sparse or not.  Returns
        SparseDataFrame if `data` is a Series or if all columns are included.
        Otherwise returns a DataFrame with some SparseBlocks.

        .. versionadded:: 0.18.2

    drop_first : bool, default False
        Whether to get k-1 dummies out of k categorical levels by removing the
        first level.

    dtype : dtype, default bool
        Data type for new columns. Only a single dtype is allowed.

        .. versionadded:: 0.18.2

    Returns
    -------
    dummies : DataFrame

    Examples
    --------
    Dask's version only works with Categorical data, as this is the only way to
    know the output shape without computing all the data.

    >>> import pandas as pd
    >>> import dask.dataframe as dd
    >>> s = dd.from_pandas(pd.Series(list('abca')), npartitions=2)
    >>> dd.get_dummies(s)
    Traceback (most recent call last):
        ...
    NotImplementedError: `get_dummies` with non-categorical dtypes is not supported...

    With categorical data:

    >>> s = dd.from_pandas(pd.Series(list('abca'), dtype='category'), npartitions=2)
    >>> dd.get_dummies(s)  # doctest: +NORMALIZE_WHITESPACE
    Dask DataFrame Structure:
                      a     b     c
    npartitions=2
    0              bool  bool  bool
    2               ...   ...   ...
    3               ...   ...   ...
    Dask Name: operation, 2 expressions
    Expr=GetDummies(frame=df)
    >>> dd.get_dummies(s).compute()  # doctest: +ELLIPSIS
           a      b      c
    0   True  False  False
    1  False   True  False
    2  False  False   True
    3   True  False  False

    See Also
    --------
    pandas.get_dummies
    """
    if isinstance(data, (pd.Series, pd.DataFrame)):
        return pd.get_dummies(
            data,
            prefix=prefix,
            prefix_sep=prefix_sep,
            dummy_na=dummy_na,
            columns=columns,
            sparse=sparse,
            drop_first=drop_first,
            dtype=dtype,
            **kwargs,
        )

    not_cat_msg = (
        "`get_dummies` with non-categorical dtypes is not "
        "supported. Please use `df.categorize()` beforehand to "
        "convert to categorical dtype."
    )

    unknown_cat_msg = (
        "`get_dummies` with unknown categories is not "
        "supported. Please use `column.cat.as_known()` or "
        "`df.categorize()` beforehand to ensure known "
        "categories"
    )

    if isinstance(data, Series):
        if not methods.is_categorical_dtype(data):
            raise NotImplementedError(not_cat_msg)
        if not has_known_categories(data):
            raise NotImplementedError(unknown_cat_msg)
    elif isinstance(data, DataFrame):
        if columns is None:
            if (data.dtypes == "object").any():
                raise NotImplementedError(not_cat_msg)
            if (data.dtypes == "string").any():
                raise NotImplementedError(not_cat_msg)
            columns = data._meta.select_dtypes(include=["category"]).columns
        else:
            if not all(methods.is_categorical_dtype(data[c]) for c in columns):
                raise NotImplementedError(not_cat_msg)

        if not all(has_known_categories(data[c]) for c in columns):
            raise NotImplementedError(unknown_cat_msg)

    return new_collection(
        GetDummies(
            data, prefix, prefix_sep, dummy_na, columns, sparse, drop_first, dtype
        )
    )


class GetDummies(Blockwise):
    _parameters = [
        "frame",
        "prefix",
        "prefix_sep",
        "dummy_na",
        "columns",
        "sparse",
        "drop_first",
        "dtype",
    ]
    _defaults = {
        "prefix": None,
        "prefix_sep": "_",
        "dummy_na": False,
        "columns": None,
        "sparse": False,
        "drop_first": False,
        "dtype": bool,
    }
    # cudf has extra kwargs after `columns`
    _keyword_only = ["sparse", "drop_first", "dtype"]

    @staticmethod
    def operation(df, *args, **kwargs):
        return get_meta_library(df).get_dummies(df, *args, **kwargs)
