from __future__ import annotations

import json
from typing import Protocol, runtime_checkable
from uuid import uuid4

import fsspec
import pandas as pd
from fsspec.implementations.local import LocalFileSystem
from packaging.version import Version

try:
    import fsspec.parquet as fsspec_parquet
except ImportError:
    fsspec_parquet = None


def _is_local_fs(fs):
    """Check if an fsspec file-system is local"""
    return fs and (
        isinstance(fs, LocalFileSystem)
        # Check wrapped pyarrow filesystem
        or _is_local_fs_pyarrow(fs)
    )


def _is_local_fs_pyarrow(fs):
    """Check if a pyarrow-based file-system is local"""
    if fs:
        if hasattr(fs, "fs"):
            # ArrowFSWrapper will have an "fs" attribute
            return _is_local_fs_pyarrow(fs.fs)
        elif hasattr(fs, "type_name"):
            # pa.fs.LocalFileSystem will have "type_name" attribute
            return fs.type_name == "local"
    return False


def _get_pyarrow_dtypes(schema, categories, dtype_backend=None):
    """Convert a pyarrow.Schema object to pandas dtype dict"""
    if dtype_backend == "numpy_nullable":
        from dask.dataframe.io.parquet.arrow import PYARROW_NULLABLE_DTYPE_MAPPING

        type_mapper = PYARROW_NULLABLE_DTYPE_MAPPING.get
    else:
        type_mapper = lambda t: t.to_pandas_dtype()

    # Check for pandas metadata
    has_pandas_metadata = schema.metadata is not None and b"pandas" in schema.metadata
    if has_pandas_metadata:
        pandas_metadata = json.loads(schema.metadata[b"pandas"].decode("utf8"))
        pandas_metadata_dtypes = {
            c.get("field_name", c.get("name", None)): c["numpy_type"]
            for c in pandas_metadata.get("columns", [])
        }
        tz = {
            c.get("field_name", c.get("name", None)): c["metadata"].get(
                "timezone", None
            )
            for c in pandas_metadata.get("columns", [])
            if c["pandas_type"] in ("datetime", "datetimetz") and c["metadata"]
        }
    else:
        pandas_metadata_dtypes = {}

    dtypes = {}
    for i in range(len(schema)):
        field = schema[i]

        # Get numpy_dtype from pandas metadata if available
        if field.name in pandas_metadata_dtypes:
            if field.name in tz:
                numpy_dtype = (
                    pd.Series([], dtype="M8[ns]").dt.tz_localize(tz[field.name]).dtype
                )
            else:
                numpy_dtype = pandas_metadata_dtypes[field.name]
        else:
            try:
                numpy_dtype = type_mapper(field.type)
            except NotImplementedError:
                continue  # Skip this field (in case we aren't reading it anyway)

        dtypes[field.name] = numpy_dtype

    if categories:
        for cat in categories:
            dtypes[cat] = "category"

    return dtypes


def _meta_from_dtypes(to_read_columns, file_dtypes, index_cols, column_index_names):
    """Get the final metadata for the dask.dataframe

    Parameters
    ----------
    to_read_columns : list
        All the columns to end up with, including index names
    file_dtypes : dict
        Mapping from column name to dtype for every element
        of ``to_read_columns``
    index_cols : list
        Subset of ``to_read_columns`` that should move to the
        index
    column_index_names : list
        The values for df.columns.name for a MultiIndex in the
        columns, or df.index.name for a regular Index in the columns

    Returns
    -------
    meta : DataFrame
    """
    data = {
        c: pd.Series([], dtype=file_dtypes.get(c, "int64")) for c in to_read_columns
    }
    indexes = [data.pop(c) for c in index_cols or []]
    if len(indexes) == 0:
        index = None
    elif len(index_cols) == 1:
        index = indexes[0]
        # XXX: this means we can't roundtrip dataframes where the index names
        # is actually __index_level_0__
        if index_cols[0] != "__index_level_0__":
            index.name = index_cols[0]
    else:
        index = pd.MultiIndex.from_arrays(indexes, names=index_cols)
    df = pd.DataFrame(data, index=index)

    if column_index_names:
        df.columns.names = column_index_names
    return df


def _guid():
    """Simple utility function to get random hex string"""
    return uuid4().hex


def _set_context(obj, stack):
    """Helper function to place an object on a context stack"""
    if stack is None:
        return obj
    return stack.enter_context(obj)


def _open_input_files(
    paths,
    fs=None,
    context_stack=None,
    open_file_func=None,
    precache_options=None,
    **kwargs,
):
    """Return a list of open-file objects given
    a list of input-file paths.

    WARNING: This utility is experimental, and is meant
    for internal ``dask.dataframe`` use only.

    Parameters
    ----------
    paths : list(str)
        Remote or local path of the parquet file
    fs : fsspec object, optional
        File-system instance to use for file handling
    context_stack : contextlib.ExitStack, Optional
        Context manager to use for open files.
    open_file_func : callable, optional
        Callable function to use for file opening. If this argument
        is specified, ``open_file_func(path, **kwargs)`` will be used
        to open each file in ``paths``. Default is ``fs.open``.
    precache_options : dict, optional
        Dictionary of key-word arguments to use for precaching.
        If ``precache_options`` contains ``{"method": "parquet"}``,
        ``fsspec.parquet.open_parquet_file`` will be used for remote
        storage.
    **kwargs :
        Key-word arguments to pass to the appropriate open function
    """
    # Use call-back function if specified
    if open_file_func is not None:
        return [
            _set_context(open_file_func(path, **kwargs), context_stack)
            for path in paths
        ]

    # Check if we are using `fsspec.parquet`.
    # In the future, fsspec should be able to handle
    # `{"method": "parquet"}`. However, for now we
    # will redirect to `open_parquet_file` manually
    precache_options = (precache_options or {}).copy()
    precache = precache_options.pop("method", None)
    if (
        precache == "parquet"
        and fs is not None
        and not _is_local_fs(fs)
        and Version(fsspec.__version__) > Version("2021.11.0")
    ):
        kwargs.update(precache_options)
        row_groups = kwargs.pop("row_groups", None) or ([None] * len(paths))
        cache_type = kwargs.pop("cache_type", "parts")
        if cache_type != "parts":
            raise ValueError(
                f"'parts' `cache_type` required for 'parquet' precaching,"
                f" got {cache_type}."
            )
        return [
            _set_context(
                fsspec_parquet.open_parquet_file(
                    path,
                    fs=fs,
                    row_groups=rgs,
                    **kwargs,
                ),
                context_stack,
            )
            for path, rgs in zip(paths, row_groups)
        ]
    elif fs is not None:
        return [_set_context(fs.open(path, **kwargs), context_stack) for path in paths]
    return [_set_context(open(path, **kwargs), context_stack) for path in paths]


@runtime_checkable
class DataFrameIOFunction(Protocol):
    """DataFrame IO function with projectable columns

    Enables column projection in ``DataFrameIOLayer``.
    """

    @property
    def columns(self):
        """Return the current column projection"""
        raise NotImplementedError

    def project_columns(self, columns):
        """Return a new DataFrameIOFunction object
        with a new column projection
        """
        raise NotImplementedError

    def __call__(self, *args, **kwargs):
        """Return a new DataFrame partition"""
        raise NotImplementedError


@runtime_checkable
class SupportsLock(Protocol):
    def acquire(self) -> object: ...

    def release(self) -> object: ...
