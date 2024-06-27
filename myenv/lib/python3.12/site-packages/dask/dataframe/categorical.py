from __future__ import annotations

from collections import defaultdict
from numbers import Integral

import pandas as pd
from pandas.api.types import is_scalar
from tlz import partition_all

from dask.base import compute_as_if_collection, tokenize
from dask.dataframe import methods
from dask.dataframe.accessor import Accessor
from dask.dataframe.dispatch import (  # noqa: F401
    categorical_dtype,
    categorical_dtype_dispatch,
    is_categorical_dtype,
)
from dask.dataframe.utils import (
    AttributeNotImplementedError,
    clear_known_categories,
    has_known_categories,
)
from dask.highlevelgraph import HighLevelGraph


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


def categorize(df, columns=None, index=None, split_every=None, **kwargs):
    """Convert columns of the DataFrame to category dtype.

    Parameters
    ----------
    columns : list, optional
        A list of column names to convert to categoricals. By default any
        column with an object dtype is converted to a categorical, and any
        unknown categoricals are made known.
    index : bool, optional
        Whether to categorize the index. By default, object indices are
        converted to categorical, and unknown categorical indices are made
        known. Set True to always categorize the index, False to never.
    split_every : int, optional
        Group partitions into groups of this size while performing a
        tree-reduction. If set to False, no tree-reduction will be used.
        Default is 16.
    kwargs
        Keyword arguments are passed on to compute.
    """
    meta = df._meta
    if columns is None:
        columns = list(meta.select_dtypes(["object", "string", "category"]).columns)
    elif is_scalar(columns):
        columns = [columns]

    # Filter out known categorical columns
    columns = [
        c
        for c in columns
        if not (is_categorical_dtype(meta[c]) and has_known_categories(meta[c]))
    ]

    if index is not False:
        if is_categorical_dtype(meta.index):
            index = not has_known_categories(meta.index)
        elif index is None:
            index = str(meta.index.dtype) in ("object", "string")

    # Nothing to do
    if not len(columns) and index is False:
        return df

    if split_every is None:
        split_every = 16
    elif split_every is False:
        split_every = df.npartitions
    elif not isinstance(split_every, Integral) or split_every < 2:
        raise ValueError("split_every must be an integer >= 2")

    token = tokenize(df, columns, index, split_every)
    a = "get-categories-chunk-" + token
    dsk = {
        (a, i): (_get_categories, key, columns, index)
        for (i, key) in enumerate(df.__dask_keys__())
    }

    prefix = "get-categories-agg-" + token
    k = df.npartitions
    depth = 0
    while k > split_every:
        b = prefix + str(depth)
        for part_i, inds in enumerate(partition_all(split_every, range(k))):
            dsk[(b, part_i)] = (_get_categories_agg, [(a, i) for i in inds])
        k = part_i + 1
        a = b
        depth += 1

    dsk[(prefix, 0)] = (_get_categories_agg, [(a, i) for i in range(k)])
    graph = HighLevelGraph.from_collections(prefix, dsk, dependencies=[df])

    # Compute the categories
    categories, index = compute_as_if_collection(
        df.__class__, graph, (prefix, 0), **kwargs
    )

    # some operations like get_dummies() rely on the order of categories
    categories = {k: v.sort_values() for k, v in categories.items()}

    # Categorize each partition
    return df.map_partitions(_categorize_block, categories, index)


class CategoricalAccessor(Accessor):
    """
    Accessor object for categorical properties of the Series values.

    Examples
    --------
    >>> s.cat.categories  # doctest: +SKIP

    Notes
    -----
    Attributes that depend only on metadata are eager

    * categories
    * ordered

    Attributes depending on the entire dataset are lazy

    * codes
    * ...

    So `df.a.cat.categories` <=> `df.a._meta.cat.categories`
    So `df.a.cat.codes` <=> `df.a.map_partitions(lambda x: x.cat.codes)`
    """

    _accessor_name = "cat"
    _accessor_methods = (
        "add_categories",
        "as_ordered",
        "as_unordered",
        "remove_categories",
        "rename_categories",
        "reorder_categories",
        "set_categories",
    )
    _accessor_properties = ()

    @property
    def known(self):
        """Whether the categories are fully known"""
        return has_known_categories(self._series)

    def as_known(self, **kwargs):
        """Ensure the categories in this series are known.

        If the categories are known, this is a no-op. If unknown, the
        categories are computed, and a new series with known categories is
        returned.

        Parameters
        ----------
        kwargs
            Keywords to pass on to the call to `compute`.
        """
        if self.known:
            return self._series
        categories = self._property_map("categories").unique().compute(**kwargs)
        return self.set_categories(categories.values)

    def as_unknown(self):
        """Ensure the categories in this series are unknown"""
        if not self.known:
            return self._series
        out = self._series.copy()
        out._meta = clear_known_categories(out._meta)
        return out

    @property
    def ordered(self):
        """Whether the categories have an ordered relationship"""
        return self._delegate_property(self._series._meta, "cat", "ordered")

    @property
    def categories(self):
        """The categories of this categorical.

        If categories are unknown, an error is raised"""
        if not self.known:
            msg = (
                "`df.column.cat.categories` with unknown categories is not "
                "supported.  Please use `column.cat.as_known()` or "
                "`df.categorize()` beforehand to ensure known categories"
            )
            raise AttributeNotImplementedError(msg)
        return self._delegate_property(self._series._meta, "cat", "categories")

    @property
    def codes(self):
        """The codes of this categorical.

        If categories are unknown, an error is raised"""
        if not self.known:
            msg = (
                "`df.column.cat.codes` with unknown categories is not "
                "supported.  Please use `column.cat.as_known()` or "
                "`df.categorize()` beforehand to ensure known categories"
            )
            raise AttributeNotImplementedError(msg)
        return self._property_map("codes")

    def remove_unused_categories(self):
        """
        Removes categories which are not used

        Notes
        -----
        This method requires a full scan of the data to compute the
        unique values, which can be expensive.
        """
        # get the set of used categories
        present = self._series.dropna().unique()
        present = pd.Index(present.compute())

        if isinstance(self._series._meta, pd.CategoricalIndex):
            meta_cat = self._series._meta
        else:
            meta_cat = self._series._meta.cat

        # Reorder to keep cat:code relationship, filtering unused (-1)
        ordered, mask = present.reindex(meta_cat.categories)
        if mask is None:
            # PANDAS-23963: old and new categories match.
            return self._series

        new_categories = ordered[mask != -1]
        meta = meta_cat.set_categories(new_categories, ordered=meta_cat.ordered)
        return self._series.map_partitions(
            self._delegate_method,
            "cat",
            "set_categories",
            (),
            {"new_categories": new_categories},
            meta=meta,
            token="cat-set_categories",
        )
