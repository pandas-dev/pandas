from __future__ import annotations

import itertools
import logging
import random
import sys
from array import array

from packaging.version import Version

from dask._compatibility import importlib_metadata
from dask.utils import Dispatch

sizeof = Dispatch(name="sizeof")

logger = logging.getLogger(__name__)


@sizeof.register(object)
def sizeof_default(o):
    return sys.getsizeof(o)


@sizeof.register(bytes)
@sizeof.register(bytearray)
def sizeof_bytes(o):
    return len(o)


@sizeof.register(memoryview)
def sizeof_memoryview(o):
    return o.nbytes


@sizeof.register(array)
def sizeof_array(o):
    return o.itemsize * len(o)


@sizeof.register(list)
@sizeof.register(tuple)
@sizeof.register(set)
@sizeof.register(frozenset)
def sizeof_python_collection(seq):
    num_items = len(seq)
    num_samples = 10
    if num_items > num_samples:
        if isinstance(seq, (set, frozenset)):
            # As of Python v3.9, it is deprecated to call random.sample() on
            # sets but since sets are unordered anyways we can simply pick
            # the first `num_samples` items.
            samples = itertools.islice(seq, num_samples)
        else:
            samples = random.sample(seq, num_samples)
        return sys.getsizeof(seq) + int(
            num_items / num_samples * sum(map(sizeof, samples))
        )
    else:
        return sys.getsizeof(seq) + sum(map(sizeof, seq))


class SimpleSizeof:
    """Sentinel class to mark a class to be skipped by the dispatcher. This only
    works if this sentinel mixin is first in the mro.

    Examples
    --------
    >>> def _get_gc_overhead():
    ...     class _CustomObject:
    ...         def __sizeof__(self):
    ...             return 0
    ...
    ...     return sys.getsizeof(_CustomObject())

    >>> class TheAnswer(SimpleSizeof):
    ...     def __sizeof__(self):
    ...         # Sizeof always add overhead of an object for GC
    ...         return 42 - _get_gc_overhead()

    >>> sizeof(TheAnswer())
    42

    """


@sizeof.register(SimpleSizeof)
def sizeof_blocked(d):
    return sys.getsizeof(d)


@sizeof.register(dict)
def sizeof_python_dict(d):
    return (
        sys.getsizeof(d)
        + sizeof(list(d.keys()))
        + sizeof(list(d.values()))
        - 2 * sizeof(list())
    )


@sizeof.register_lazy("cupy")
def register_cupy():
    import cupy

    @sizeof.register(cupy.ndarray)
    def sizeof_cupy_ndarray(x):
        return int(x.nbytes)


@sizeof.register_lazy("numba")
def register_numba():
    import numba.cuda

    @sizeof.register(numba.cuda.cudadrv.devicearray.DeviceNDArray)
    def sizeof_numba_devicendarray(x):
        return int(x.nbytes)


@sizeof.register_lazy("rmm")
def register_rmm():
    import rmm

    # Only included in 0.11.0+
    if hasattr(rmm, "DeviceBuffer"):

        @sizeof.register(rmm.DeviceBuffer)
        def sizeof_rmm_devicebuffer(x):
            return int(x.nbytes)


@sizeof.register_lazy("numpy")
def register_numpy():
    import numpy as np

    @sizeof.register(np.ndarray)
    def sizeof_numpy_ndarray(x):
        if 0 in x.strides:
            xs = x[tuple(slice(None) if s != 0 else slice(1) for s in x.strides)]
            return xs.nbytes
        return int(x.nbytes)


@sizeof.register_lazy("pandas")
def register_pandas():
    import numpy as np
    import pandas as pd

    OBJECT_DTYPES = (object, pd.StringDtype("python"))

    def object_size(*xs):
        if not xs:
            return 0
        ncells = sum(len(x) for x in xs)
        if not ncells:
            return 0

        # Deduplicate Series of references to the same objects,
        # e.g. as produced by read_parquet
        unique_samples = {}
        for x in xs:
            sample = np.random.choice(x, size=100, replace=True)
            for i in sample.tolist():
                unique_samples[id(i)] = i

        nsamples = 100 * len(xs)
        sample_nbytes = sum(sizeof(i) for i in unique_samples.values())
        if len(unique_samples) / nsamples > 0.5:
            # Less than half of the references are duplicated.
            # Assume that, if we were to analyze twice the amount of random references,
            # we would get twice the amount of unique objects too.
            return int(sample_nbytes * ncells / nsamples)
        else:
            # Assume we've already found all unique objects and that all references that
            # we have not yet analyzed are going to point to the same data.
            return sample_nbytes

    @sizeof.register(pd.DataFrame)
    def sizeof_pandas_dataframe(df):
        p = sizeof(df.index) + sizeof(df.columns)
        object_cols = []
        prev_dtype = None

        # Unlike df.items(), df._series will not duplicate multiple views of the same
        # column e.g. df[["x", "x", "x"]]
        for col in df._series.values():
            if prev_dtype is None or col.dtype != prev_dtype:
                prev_dtype = col.dtype
                # Contiguous columns of the same dtype share the same overhead
                p += 1200
            p += col.memory_usage(index=False, deep=False)
            if col.dtype in OBJECT_DTYPES:
                object_cols.append(col._values)

        # Deduplicate references to the same objects appearing in different Series
        p += object_size(*object_cols)

        return max(1200, p)

    @sizeof.register(pd.Series)
    def sizeof_pandas_series(s):
        # https://github.com/dask/dask/pull/9776#issuecomment-1359085962
        p = 1200 + sizeof(s.index) + s.memory_usage(index=False, deep=False)
        if s.dtype in OBJECT_DTYPES:
            p += object_size(s._values)
        return p

    @sizeof.register(pd.Index)
    def sizeof_pandas_index(i):
        p = 400 + i.memory_usage(deep=False)
        if i.dtype in OBJECT_DTYPES:
            p += object_size(i)
        return p

    @sizeof.register(pd.MultiIndex)
    def sizeof_pandas_multiindex(i):
        return sum(sizeof(l) for l in i.levels) + sum(c.nbytes for c in i.codes)


@sizeof.register_lazy("scipy")
def register_spmatrix():
    import scipy
    from scipy import sparse

    if Version(scipy.__version__) < Version("1.12.0.dev0"):

        @sizeof.register(sparse.dok_matrix)
        def sizeof_spmatrix_dok(s):
            return s.__sizeof__()

    @sizeof.register(sparse.spmatrix)
    def sizeof_spmatrix(s):
        return sum(sizeof(v) for v in s.__dict__.values())


@sizeof.register_lazy("pyarrow")
def register_pyarrow():
    import pyarrow as pa

    def _get_col_size(data):
        p = 0
        if not isinstance(data, pa.ChunkedArray):
            data = data.data  # pyarrow <0.15.0
        for chunk in data.iterchunks():
            for buffer in chunk.buffers():
                if buffer:
                    p += buffer.size
        return p

    @sizeof.register(pa.Table)
    def sizeof_pyarrow_table(table):
        p = sizeof(table.schema.metadata)
        for col in table.itercolumns():
            p += _get_col_size(col)
        return int(p) + 1000

    @sizeof.register(pa.ChunkedArray)
    def sizeof_pyarrow_chunked_array(data):
        return int(_get_col_size(data)) + 1000


@sizeof.register_lazy("xarray")
def register_xarray():
    import sys

    import xarray as xr

    XARRAY_VERSION = Version(xr.__version__)
    XARRAY_GE_2024_02 = XARRAY_VERSION >= Version("2024.02.0")

    @sizeof.register(xr.core.utils.Frozen)
    def xarray_sizeof_frozen(obj):
        return sys.getsizeof(obj) + sizeof(obj.mapping)

    @sizeof.register(xr.DataArray)
    @sizeof.register(xr.Variable)
    def xarray_sizeof_da(obj):
        return sys.getsizeof(obj) + sizeof(obj.data)

    @sizeof.register(xr.Dataset)
    def xarray_sizeof_ds(obj):
        return sys.getsizeof(obj) + sizeof(obj.variables)

    if XARRAY_GE_2024_02:
        xarray_sizeof_da = sizeof.register(xr.NamedArray)(xarray_sizeof_da)

    @sizeof.register(xr.core.indexes.Indexes)
    def xarray_sizeof_indexes(obj):
        return (
            sys.getsizeof(obj)
            + sizeof(obj._index_type)
            + sizeof(obj._indexes)
            + sizeof(obj._variables)
            + sizeof(obj._dims)
        )

    @sizeof.register(xr.core.indexes.PandasIndex)
    def xarray_sizeof_pd_index(obj):
        return (
            sys.getsizeof(obj)
            + sizeof(obj.index)
            + sizeof(obj.dim)
            + sizeof(obj.coord_dtype)
        )


def _register_entry_point_plugins():
    """Register sizeof implementations exposed by the entry_point mechanism."""
    for entry_point in importlib_metadata.entry_points(group="dask.sizeof"):
        registrar = entry_point.load()
        try:
            registrar(sizeof)
        except Exception:
            logger.exception(
                f"Failed to register sizeof entry point {entry_point.name}"
            )


_register_entry_point_plugins()
