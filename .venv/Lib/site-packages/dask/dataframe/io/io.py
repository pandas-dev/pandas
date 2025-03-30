from __future__ import annotations

import bisect
from threading import Lock
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from dask.blockwise import BlockwiseDepDict, blockwise
from dask.dataframe.dispatch import meta_lib_from_array, tolist
from dask.dataframe.utils import pyarrow_strings_enabled
from dask.highlevelgraph import HighLevelGraph
from dask.tokenize import tokenize
from dask.utils import ensure_dict, is_series_like

if TYPE_CHECKING:
    pass


lock = Lock()


def _meta_from_array(x, columns=None, index=None, meta=None):
    """Create empty DataFrame or Series which has correct dtype"""

    if x.ndim > 2:
        raise ValueError(
            "from_array does not input more than 2D array, got"
            " array with shape %r" % (x.shape,)
        )

    if index is not None:
        from dask.dataframe import Index

        if not isinstance(index, Index):
            raise ValueError("'index' must be an instance of dask.dataframe.Index")
        index = index._meta

    if meta is None:
        meta = meta_lib_from_array(x).DataFrame()

    if getattr(x.dtype, "names", None) is not None:
        # record array has named columns
        if columns is None:
            columns = list(x.dtype.names)
        elif np.isscalar(columns):
            raise ValueError("For a struct dtype, columns must be a list.")
        elif not all(i in x.dtype.names for i in columns):
            extra = sorted(set(columns).difference(x.dtype.names))
            raise ValueError(f"dtype {x.dtype} doesn't have fields {extra}")
        fields = x.dtype.fields
        dtypes = [fields[n][0] if n in fields else "f8" for n in columns]
    elif x.ndim == 1:
        if np.isscalar(columns) or columns is None:
            return meta._constructor_sliced(
                [], name=columns, dtype=x.dtype, index=index
            )
        elif len(columns) == 1:
            return meta._constructor(
                np.array([], dtype=x.dtype), columns=columns, index=index
            )
        raise ValueError(
            "For a 1d array, columns must be a scalar or single element list"
        )
    else:
        if np.isnan(x.shape[1]):
            raise ValueError("Shape along axis 1 must be known")
        if columns is None:
            columns = list(range(x.shape[1])) if x.ndim == 2 else [0]
        elif len(columns) != x.shape[1]:
            raise ValueError(
                "Number of column names must match width of the array. "
                f"Got {len(columns)} names for {x.shape[1]} columns"
            )
        dtypes = [x.dtype] * len(columns)

    data = {c: np.array([], dtype=dt) for (c, dt) in zip(columns, dtypes)}
    return meta._constructor(data, columns=columns, index=index)


def _partition_from_array(data, index=None, initializer=None, **kwargs):
    """Create a Dask partition for either a DataFrame or Series.

    Designed to be used with :func:`dask.blockwise.blockwise`. ``data`` is the array
    from which the partition will be created. ``index`` can be:

    1. ``None``, in which case each partition has an independent RangeIndex
    2. a `tuple` with two elements, the start and stop values for a RangeIndex for
       this partition, which gives a continuously varying RangeIndex over the
       whole Dask DataFrame
    3. an instance of a ``pandas.Index`` or a subclass thereof

    The ``kwargs`` _must_ contain an ``initializer`` key which is set by calling
    ``type(meta)``.
    """
    if isinstance(index, tuple):
        index = pd.RangeIndex(*index)
    return initializer(data, index=index, **kwargs)


def from_dask_array(x, columns=None, index=None, meta=None):
    """Create a Dask DataFrame from a Dask Array.

    Converts a 2d array into a DataFrame and a 1d array into a Series.

    Parameters
    ----------
    x : da.Array
    columns : list or string
        list of column names if DataFrame, single string if Series
    index : dask.dataframe.Index, optional
        An optional *dask* Index to use for the output Series or DataFrame.

        The default output index depends on whether `x` has any unknown
        chunks. If there are any unknown chunks, the output has ``None``
        for all the divisions (one per chunk). If all the chunks are known,
        a default index with known divisions is created.

        Specifying `index` can be useful if you're conforming a Dask Array
        to an existing dask Series or DataFrame, and you would like the
        indices to match.
    meta : object, optional
        An optional `meta` parameter can be passed for dask
        to specify the concrete dataframe type to be returned.
        By default, pandas DataFrame is used.

    Examples
    --------
    >>> import dask.array as da
    >>> import dask.dataframe as dd
    >>> x = da.ones((4, 2), chunks=(2, 2))
    >>> df = dd.io.from_dask_array(x, columns=['a', 'b'])
    >>> df.compute()
         a    b
    0  1.0  1.0
    1  1.0  1.0
    2  1.0  1.0
    3  1.0  1.0

    See Also
    --------
    dask.bag.to_dataframe: from dask.bag
    """
    meta = _meta_from_array(x, columns, index, meta=meta)

    name = "from-dask-array-" + tokenize(x, columns)
    graph_dependencies = [x]
    arrays_and_indices = [x.name, "ij" if x.ndim == 2 else "i"]
    numblocks = {x.name: x.numblocks}

    if index is not None:
        # An index is explicitly given by the caller, so we can pass it through to the
        # initializer after a few checks.
        if index.npartitions != x.numblocks[0]:
            msg = (
                "The index and array have different numbers of blocks. "
                "({} != {})".format(index.npartitions, x.numblocks[0])
            )
            raise ValueError(msg)
        divisions = index.divisions
        graph_dependencies.append(index)
        arrays_and_indices.extend([index._name, "i"])
        numblocks[index._name] = (index.npartitions,)

    elif np.isnan(sum(x.shape)):
        # The shape of the incoming array is not known in at least one dimension. As
        # such, we can't create an index for the entire output DataFrame and we set
        # the divisions to None to represent that.
        divisions = [None] * (len(x.chunks[0]) + 1)
    else:
        # The shape of the incoming array is known and we don't have an explicit index.
        # Create a mapping of chunk number in the incoming array to
        # (start row, stop row) tuples. These tuples will be used to create a sequential
        # RangeIndex later on that is continuous over the whole DataFrame.
        n_elements = sum(x.chunks[0])
        divisions = [0]
        stop = 0
        index_mapping = {}
        for i, increment in enumerate(x.chunks[0]):
            stop += increment
            index_mapping[(i,)] = (divisions[i], stop)

            # last division corrected, even if there are empty chunk(s) at the end
            if stop == n_elements:
                stop -= 1

            divisions.append(stop)

        arrays_and_indices.extend([BlockwiseDepDict(mapping=index_mapping), "i"])

    if is_series_like(meta):
        kwargs = {"dtype": x.dtype, "name": meta.name, "initializer": type(meta)}
    else:
        kwargs = {"columns": meta.columns, "initializer": type(meta)}

    blk = blockwise(
        _partition_from_array,
        name,
        "i",
        *arrays_and_indices,
        numblocks=numblocks,
        concatenate=True,
        # kwargs passed through to the DataFrame/Series initializer
        **kwargs,
    )

    graph = HighLevelGraph.from_collections(name, blk, dependencies=graph_dependencies)

    from dask.array.optimization import optimize
    from dask.dataframe.dask_expr._collection import from_graph, new_collection
    from dask.dataframe.dask_expr._expr import ArrowStringConversion
    from dask.utils import key_split

    keys = [(name, i) for i in range(len(divisions) - 1)]
    result = from_graph(
        optimize(ensure_dict(graph), keys),
        meta,
        divisions,
        keys,
        key_split(name),
    )
    if pyarrow_strings_enabled():
        return new_collection(ArrowStringConversion(result.expr))
    return result


def _link(token, result):
    """A dummy function to link results together in a graph

    We use this to enforce an artificial sequential ordering on tasks that
    don't explicitly pass around a shared resource
    """
    return None


def _df_to_bag(df, index=False, format="tuple"):
    if isinstance(df, pd.DataFrame):
        if format == "tuple":
            return list(map(tuple, df.itertuples(index)))
        elif format == "dict":
            if index:
                return [
                    {**{"index": idx}, **values}
                    for values, idx in zip(df.to_dict("records"), df.index)
                ]
            else:
                return df.to_dict(orient="records")
    elif isinstance(df, pd.Series):
        if format == "tuple":
            return list(df.items()) if index else list(df)
        elif format == "dict":
            return df.to_frame().to_dict(orient="records")


def sorted_division_locations(seq, npartitions=None, chunksize=None):
    """Find division locations and values in sorted list

    Examples
    --------

    >>> L = ['A', 'B', 'C', 'D', 'E', 'F']
    >>> sorted_division_locations(L, chunksize=2)
    (['A', 'C', 'E', 'F'], [0, 2, 4, 6])

    >>> sorted_division_locations(L, chunksize=3)
    (['A', 'D', 'F'], [0, 3, 6])

    >>> L = ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'C']
    >>> sorted_division_locations(L, chunksize=3)
    (['A', 'B', 'C', 'C'], [0, 4, 7, 8])

    >>> sorted_division_locations(L, chunksize=2)
    (['A', 'B', 'C', 'C'], [0, 4, 7, 8])

    >>> sorted_division_locations(['A'], chunksize=2)
    (['A', 'A'], [0, 1])
    """
    if (npartitions is None) == (chunksize is None):
        raise ValueError("Exactly one of npartitions and chunksize must be specified.")

    # Convert from an ndarray to a plain list so that
    # any divisions we extract from seq are plain Python scalars.
    seq = tolist(seq)
    # we use bisect later, so we need sorted.
    seq_unique = sorted(set(seq))
    duplicates = len(seq_unique) < len(seq)
    enforce_exact = False

    if duplicates:
        offsets = [bisect.bisect_left(seq, x) for x in seq_unique]
        enforce_exact = npartitions and len(offsets) >= npartitions
    else:
        offsets = seq_unique = None

    # Define chunksize and residual so that
    # npartitions can be exactly satisfied
    # when duplicates is False
    residual = 0
    subtract_drift = False
    if npartitions:
        chunksize = len(seq) // npartitions
        residual = len(seq) % npartitions
        subtract_drift = True

    def chunksizes(ind):
        # Helper function to satisfy npartitions
        return chunksize + int(ind < residual)

    # Always start with 0th item in seqarr,
    # and then try to take chunksize steps
    # along the seqarr array
    divisions = [seq[0]]
    locations = [0]
    i = chunksizes(0)
    ind = None  # ind cache (sometimes avoids nonzero call)
    drift = 0  # accumulated drift away from ideal chunksizes
    divs_remain = npartitions - len(divisions) if enforce_exact else None
    while i < len(seq):
        # Map current position selection (i)
        # to the corresponding division value (div)
        div = seq[i]
        # pos is the position of the first occurrence of
        # div (which is i when seq has no duplicates)
        if duplicates:
            # Note: cupy requires casts to `int` below
            if ind is None:
                ind = bisect.bisect_left(seq_unique, seq[i])
            if enforce_exact:
                # Avoid "over-stepping" too many unique
                # values when npartitions is approximately
                # equal to len(offsets)
                offs_remain = len(offsets) - ind
                if divs_remain > offs_remain:
                    ind -= divs_remain - offs_remain
                    i = offsets[ind]
                    div = seq[i]
            pos = int(offsets[ind])
        else:
            pos = i
        if div <= divisions[-1]:
            # pos overlaps with divisions.
            # Try the next element on the following pass
            if duplicates:
                ind += 1
                # Note: cupy requires cast to `int`
                i = int(offsets[ind]) if ind < len(offsets) else len(seq)
            else:
                i += 1
        else:
            # pos does not overlap with divisions.
            # Append candidate pos/div combination, and
            # take another chunksize step
            if subtract_drift:
                # Only subtract drift when user specified npartitions
                drift = drift + ((pos - locations[-1]) - chunksizes(len(divisions) - 1))
            if enforce_exact:
                divs_remain -= 1
            i = pos + max(1, chunksizes(len(divisions)) - drift)
            divisions.append(div)
            locations.append(pos)
            ind = None

    # The final element of divisions/locations
    # will always be the same
    divisions.append(seq[-1])
    locations.append(len(seq))

    return divisions, locations
