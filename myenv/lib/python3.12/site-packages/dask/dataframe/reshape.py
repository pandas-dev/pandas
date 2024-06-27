from __future__ import annotations

import numpy as np
import pandas as pd
from pandas.api.types import is_list_like, is_scalar

import dask
from dask.dataframe import methods
from dask.dataframe._compat import PANDAS_GE_200
from dask.dataframe.core import DataFrame, Series, apply_concat_apply, map_partitions
from dask.dataframe.utils import has_known_categories
from dask.typing import no_default
from dask.utils import M, get_meta_library

###############################################################
# Dummies
###############################################################


_get_dummies_dtype_default = bool if PANDAS_GE_200 else np.uint8


def get_dummies(
    data,
    prefix=None,
    prefix_sep="_",
    dummy_na=False,
    columns=None,
    sparse=False,
    drop_first=False,
    dtype=_get_dummies_dtype_default,
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
                       a      b      c
    npartitions=2
    0              bool  bool  bool
    2                ...    ...    ...
    3                ...    ...    ...
    Dask Name: get_dummies, 2 graph layers
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

    return map_partitions(
        get_meta_library(data).get_dummies,
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


###############################################################
# Pivot table
###############################################################


def pivot_table(df, index=None, columns=None, values=None, aggfunc="mean"):
    """
    Create a spreadsheet-style pivot table as a DataFrame. Target ``columns``
    must have category dtype to infer result's ``columns``.
    ``index``, ``columns``, and ``aggfunc`` must be all scalar.
    ``values`` can be scalar or list-like.

    Parameters
    ----------
    df : DataFrame
    index : scalar
        column to be index
    columns : scalar
        column to be columns
    values : scalar or list(scalar)
        column(s) to aggregate
    aggfunc : {'mean', 'sum', 'count', 'first', 'last'}, default 'mean'

    Returns
    -------
    table : DataFrame

    See Also
    --------
    pandas.DataFrame.pivot_table
    """

    if not is_scalar(index) or index is None:
        raise ValueError("'index' must be the name of an existing column")
    if not is_scalar(columns) or columns is None:
        raise ValueError("'columns' must be the name of an existing column")
    if not methods.is_categorical_dtype(df[columns]):
        raise ValueError("'columns' must be category dtype")
    if not has_known_categories(df[columns]):
        raise ValueError(
            "'columns' must have known categories. Please use "
            "`df[columns].cat.as_known()` beforehand to ensure "
            "known categories"
        )
    if not (
        is_list_like(values)
        and all([is_scalar(v) for v in values])
        or is_scalar(values)
    ):
        raise ValueError("'values' must refer to an existing column or columns")

    available_aggfuncs = ["mean", "sum", "count", "first", "last"]

    if not is_scalar(aggfunc) or aggfunc not in available_aggfuncs:
        raise ValueError(
            "aggfunc must be either " + ", ".join(f"'{x}'" for x in available_aggfuncs)
        )

    # _emulate can't work for empty data
    # the result must have CategoricalIndex columns

    columns_contents = pd.CategoricalIndex(df[columns].cat.categories, name=columns)
    if is_scalar(values):
        new_columns = columns_contents
    else:
        new_columns = pd.MultiIndex.from_product(
            (sorted(values), columns_contents), names=[None, columns]
        )

    if aggfunc in ["first", "last"]:
        # Infer datatype as non-numeric values are allowed
        if is_scalar(values):
            meta = pd.DataFrame(
                columns=new_columns,
                dtype=df[values].dtype,
                index=pd.Index(df._meta[index]),
            )
        else:
            meta = pd.DataFrame(
                columns=new_columns,
                index=pd.Index(df._meta[index]),
            )
            for value_col in values:
                meta[value_col] = meta[value_col].astype(df[values].dtypes[value_col])
    else:
        # Use float64 as other aggregate functions require numerical data
        meta = pd.DataFrame(
            columns=new_columns, dtype=np.float64, index=pd.Index(df._meta[index])
        )

    kwargs = {"index": index, "columns": columns, "values": values}

    if aggfunc in ["sum", "mean"]:
        pv_sum = apply_concat_apply(
            [df],
            chunk=methods.pivot_sum,
            aggregate=methods.pivot_agg,
            meta=meta,
            token="pivot_table_sum",
            chunk_kwargs=kwargs,
        )

    if aggfunc in ["count", "mean"]:
        pv_count = apply_concat_apply(
            [df],
            chunk=methods.pivot_count,
            aggregate=methods.pivot_agg,
            meta=meta,
            token="pivot_table_count",
            chunk_kwargs=kwargs,
        )

    if aggfunc == "sum":
        return pv_sum
    elif aggfunc == "count":
        return pv_count
    elif aggfunc == "mean":
        return pv_sum / pv_count
    elif aggfunc == "first":
        return apply_concat_apply(
            [df],
            chunk=methods.pivot_first,
            aggregate=methods.pivot_agg_first,
            meta=meta,
            token="pivot_table_first",
            chunk_kwargs=kwargs,
        )
    elif aggfunc == "last":
        return apply_concat_apply(
            [df],
            chunk=methods.pivot_last,
            aggregate=methods.pivot_agg_last,
            meta=meta,
            token="pivot_table_last",
            chunk_kwargs=kwargs,
        )
    else:
        raise ValueError


###############################################################
# Melt
###############################################################


def melt(
    frame,
    id_vars=None,
    value_vars=None,
    var_name=None,
    value_name="value",
    col_level=None,
):
    """
    Unpivots a DataFrame from wide format to long format, optionally leaving identifier variables set.

    This function is useful to massage a DataFrame into a format where one or more columns are identifier variables
    (``id_vars``), while all other columns, considered measured variables (``value_vars``), are "unpivoted" to the row
    axis, leaving just two non-identifier columns, 'variable' and 'value'.

    Parameters
    ----------
    frame : DataFrame
    id_vars : tuple, list, or ndarray, optional
        Column(s) to use as identifier variables.
    value_vars : tuple, list, or ndarray, optional
        Column(s) to unpivot. If not specified, uses all columns that
        are not set as `id_vars`.
    var_name : scalar
        Name to use for the 'variable' column. If None it uses
        ``frame.columns.name`` or 'variable'.
    value_name : scalar, default 'value'
        Name to use for the 'value' column.
    col_level : int or string, optional
        If columns are a MultiIndex then use this level to melt.

    Returns
    -------
    DataFrame
        Unpivoted DataFrame.

    See Also
    --------
    pandas.DataFrame.melt
    """
    # let pandas do upcasting as needed during melt
    with dask.config.set({"dataframe.convert-string": False}):
        return frame.map_partitions(
            M.melt,
            meta=no_default,
            id_vars=id_vars,
            value_vars=value_vars,
            var_name=var_name,
            value_name=value_name,
            col_level=col_level,
            token="melt",
        )
