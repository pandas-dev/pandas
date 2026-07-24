from __future__ import annotations

from collections import defaultdict

from dask.dataframe import methods
from dask.dataframe.dispatch import (  # noqa: F401
    categorical_dtype,
    categorical_dtype_dispatch,
    is_categorical_dtype,
)


def _categorize_block(df, categories, index):
    """Categorize a dataframe with given categories

    df: DataFrame
    categories: dict mapping column name to iterable of categories
    """
    df = df.copy()
    for col, vals in categories.items():
        if is_categorical_dtype(df[col]):
            df[col] = df[col].cat.set_categories(vals)
        else:
            cat_dtype = categorical_dtype(meta=df[col], categories=vals, ordered=False)
            df[col] = df[col].astype(cat_dtype)
    if index is not None:
        if is_categorical_dtype(df.index):
            ind = df.index.set_categories(index)
        else:
            cat_dtype = categorical_dtype(
                meta=df.index, categories=index, ordered=False
            )
            ind = df.index.astype(dtype=cat_dtype)
        ind.name = df.index.name
        df.index = ind
    return df


def _get_categories(df, columns, index):
    res = {}
    for col in columns:
        x = df[col]
        if is_categorical_dtype(x):
            res[col] = x._constructor(x.cat.categories)
        else:
            res[col] = x.dropna().drop_duplicates()
    if index:
        if is_categorical_dtype(df.index):
            return res, df.index.categories
        return res, df.index.dropna().drop_duplicates()
    return res, None


def _get_categories_agg(parts):
    res = defaultdict(list)
    res_ind = []
    for p in parts:
        for k, v in p[0].items():
            res[k].append(v)
        res_ind.append(p[1])
    res = {
        k: methods.concat(v, ignore_index=True).drop_duplicates()
        for k, v in res.items()
    }
    if res_ind[0] is None:
        return res, None
    return res, res_ind[0].append(res_ind[1:]).drop_duplicates()
