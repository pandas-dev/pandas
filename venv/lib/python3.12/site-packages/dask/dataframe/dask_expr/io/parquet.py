from __future__ import annotations

import contextlib
import itertools
import operator
import os
import pickle
import statistics
import warnings
import weakref
from abc import abstractmethod
from collections import defaultdict
from functools import cached_property, partial

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset
import pyarrow.dataset as pa_ds
import pyarrow.fs as pa_fs
import pyarrow.parquet as pq
import tlz as toolz
from fsspec.utils import stringify_path
from toolz import identity

import dask
from dask._task_spec import Task, TaskRef
from dask.core import flatten
from dask.dataframe.dask_expr._expr import (
    EQ,
    GE,
    GT,
    LE,
    LT,
    NE,
    And,
    Blockwise,
    Expr,
    Filter,
    Index,
    Lengths,
    Literal,
    Or,
    Projection,
    determine_column_projection,
)
from dask.dataframe.dask_expr._reductions import Len
from dask.dataframe.dask_expr._util import _convert_to_list
from dask.dataframe.dask_expr.io import BlockwiseIO, PartitionsFiltered
from dask.dataframe.dask_expr.io.io import FusedParquetIO
from dask.dataframe.io.parquet.core import (
    ParquetFunctionWrapper,
    ToParquetFunctionWrapper,
    aggregate_row_groups,
    apply_filters,
    get_engine,
    set_index_columns,
    sorted_columns,
)
from dask.dataframe.io.parquet.utils import _split_user_options
from dask.dataframe.io.utils import _is_local_fs
from dask.delayed import delayed
from dask.tokenize import _tokenize_deterministic, normalize_token, tokenize
from dask.typing import Key
from dask.utils import apply, funcname, natural_sort_key, parse_bytes, typename


@normalize_token.register(pa.fs.FileInfo)
def _tokenize_fileinfo(fileinfo):
    return type(fileinfo).__name__, (
        fileinfo.path,
        fileinfo.size,
        fileinfo.mtime_ns,
    )


_CPU_COUNT_SET = False


def _maybe_adjust_cpu_count():
    global _CPU_COUNT_SET
    if not _CPU_COUNT_SET:
        # Set the number of threads to the number of cores
        # This is a default for pyarrow, but it's not set by default in
        # dask/distributed
        pa.set_cpu_count(os.cpu_count())
        _CPU_COUNT_SET = True


_STATS_CACHE = {}  # type: ignore


PYARROW_NULLABLE_DTYPE_MAPPING = {
    pa.int8(): pd.Int8Dtype(),
    pa.int16(): pd.Int16Dtype(),
    pa.int32(): pd.Int32Dtype(),
    pa.int64(): pd.Int64Dtype(),
    pa.uint8(): pd.UInt8Dtype(),
    pa.uint16(): pd.UInt16Dtype(),
    pa.uint32(): pd.UInt32Dtype(),
    pa.uint64(): pd.UInt64Dtype(),
    pa.bool_(): pd.BooleanDtype(),
    pa.string(): pd.StringDtype(),
    pa.float32(): pd.Float32Dtype(),
    pa.float64(): pd.Float64Dtype(),
}

NONE_LABEL = "__null_dask_index__"

_CACHED_PLAN_SIZE = 10
_cached_plan = {}  # type: ignore


class FragmentWrapper:
    _filesystems = weakref.WeakValueDictionary()  # type: ignore
    _filesystem_pickle_cache = (-1, None)

    def __init__(
        self, fragment, filesystem, file_size=None, fragment_packed=None
    ) -> None:
        """Wrap a pyarrow Fragment to only deserialize when needed."""
        # https://github.com/apache/arrow/issues/40279
        self._fragment = fragment
        self._fragment_packed = fragment_packed
        self._file_size = file_size
        self._fs = None
        self._filesystem = filesystem

    def pack(self):
        if self._fragment_packed is None:
            part_expr = self._fragment.partition_expression
            if part_expr.equals(pc.scalar(True)):
                part_expr = True
            pqformat = self._fragment.format
            if pqformat.read_options.equals(pa_ds.ParquetFileFormat().read_options):
                pqformat = None

            fs = self._filesystem or self._fragment.filesystem
            assert fs.equals(self._fragment.filesystem)
            if self._filesystem_pickle_cache[0] != id(fs):
                FragmentWrapper._filesystem_pickle_cache = (id(fs), pickle.dumps(fs))
            fs_pkl = FragmentWrapper._filesystem_pickle_cache[1]
            self._fragment_packed = (
                pqformat,
                (
                    self._fragment.path
                    if self._fragment.buffer is None
                    else self._fragment.buffer
                ),
                fs_pkl,
                part_expr,
                self._file_size,
            )
        self._fs = self._fragment = None

    def unpack(self):
        if self._fragment is None:
            (
                pqformat,
                path_or_buffer,
                fs_raw,
                partition_expression,
                file_size,
            ) = self._fragment_packed
            fs = FragmentWrapper._filesystems.get(fs_raw)
            if fs is None:
                fs = pickle.loads(fs_raw)
                FragmentWrapper._filesystems[fs_raw] = fs
            if partition_expression is True:
                partition_expression = pc.scalar(True)
            # arrow doesn't keep the python object alive so if we want to reuse
            # we need to keep a reference
            self._fs = fs
            if pqformat is None:
                pqformat = pa_ds.ParquetFileFormat()
            self._fragment = pqformat.make_fragment(
                path_or_buffer,
                filesystem=fs,
                partition_expression=partition_expression,
                file_size=file_size,
            )
        self._fragment_packed = None

    @property
    def fragment(self):
        self.unpack()
        return self._fragment

    def __dask_tokenize__(self):
        return type(self).__name__, normalize_token(
            (
                self.fragment,
                self._fragment_packed,
                self._file_size,
            )
        )

    def __reduce__(self):
        self.pack()
        return FragmentWrapper, (None, None, None, self._fragment_packed)


def _control_cached_plan(key):
    if len(_cached_plan) > _CACHED_PLAN_SIZE and key not in _cached_plan:
        key_to_pop = list(_cached_plan.keys())[0]
        _cached_plan.pop(key_to_pop)


@normalize_token.register(pa_ds.Dataset)
def normalize_pa_ds(ds):
    return (ds.files, ds.schema)


@normalize_token.register(pa_ds.FileFormat)
def normalize_pa_file_format(file_format):
    return str(file_format)


@normalize_token.register(pa.Schema)
def normalize_pa_schema(schema):
    return schema.to_string()


@normalize_token.register(pq.ParquetSchema)
def normalize_pq_schema(schema):
    try:
        return hash(schema)
    except TypeError:  # pyarrow version not supporting ParquetSchema hash
        return hash(repr(schema))


@normalize_token.register(pq.FileMetaData)
def normalize_pq_filemetadata(meta):
    try:
        return hash(meta)
    except TypeError:
        # pyarrow version not implementing hash for FileMetaData
        # use same logic as implemented in version that does support hashing
        # https://github.com/apache/arrow/blob/bbe59b35de33a0534fc76c9617aa4746031ce16c/python/pyarrow/_parquet.pyx#L853
        return hash(
            (
                repr(meta.schema),
                meta.num_rows,
                meta.num_row_groups,
                meta.format_version,
                meta.serialized_size,
            )
        )


class ToParquet(Expr):
    _parameters = [
        "frame",
        "path",
        "fs",
        "fmd",
        "engine",
        "offset",
        "partition_on",
        "write_metadata_file",
        "name_function",
        "write_kwargs",
        "append",
    ]

    @property
    def _meta(self):
        return None

    def _divisions(self):
        return (None, None)

    def _lower(self):
        return ToParquetBarrier(
            ToParquetData(
                *self.operands,
            ),
            *self.operands[1:],
        )


class ToParquetData(Blockwise):
    _parameters = ToParquet._parameters

    @property
    def io_func(self):
        return ToParquetFunctionWrapper(
            self.engine,
            self.path,
            self.fs,
            self.partition_on,
            self.write_metadata_file,
            self.offset,
            self.name_function,
            self.write_kwargs,
        )

    def _divisions(self):
        return (None,) * (self.frame.npartitions + 1)

    def _task(self, name: Key, index: int) -> Task:
        return Task(name, self.io_func, TaskRef((self.frame._name, index)), (index,))


class ToParquetBarrier(Expr):
    _parameters = ToParquet._parameters

    @property
    def _meta(self):
        return None

    def _divisions(self):
        return (None, None)

    def _layer(self):
        if self.write_metadata_file:
            append = self.append
            compression = self.write_kwargs.get("compression")
            return {
                (self._name, 0): (
                    apply,
                    self.engine.write_metadata,
                    [
                        self.frame.__dask_keys__(),
                        self.fmd,
                        self.fs,
                        self.path,
                    ],
                    {"append": append, "compression": compression},
                )
            }
        else:
            return {(self._name, 0): (lambda x: None, self.frame.__dask_keys__())}


def to_parquet(
    df,
    path,
    compression="snappy",
    write_index=True,
    append=False,
    overwrite=False,
    ignore_divisions=False,
    partition_on=None,
    storage_options=None,
    custom_metadata=None,
    write_metadata_file=None,
    compute=True,
    compute_kwargs=None,
    schema="infer",
    name_function=None,
    filesystem=None,
    engine=None,
    **kwargs,
):
    """Store Dask.dataframe to Parquet files

    Notes
    -----
    Each partition will be written to a separate file.

    Parameters
    ----------
    df : dask.dataframe.DataFrame
    path : string or pathlib.Path
        Destination directory for data.  Prepend with protocol like ``s3://``
        or ``hdfs://`` for remote data.
    compression : string or dict, default 'snappy'
        Either a string like ``"snappy"`` or a dictionary mapping column names
        to compressors like ``{"name": "gzip", "values": "snappy"}``. Defaults
        to ``"snappy"``.
    write_index : bool, default True
        Whether or not to write the index. Defaults to True.
    append : bool, default False
        If False (default), construct data-set from scratch. If True, add new
        row-group(s) to an existing data-set. In the latter case, the data-set
        must exist, and the schema must match the input data.
    overwrite : bool, default False
        Whether or not to remove the contents of `path` before writing the dataset.
        The default is False.  If True, the specified path must correspond to
        a directory (but not the current working directory).  This option cannot
        be set to True if `append=True`.
        NOTE: `overwrite=True` will remove the original data even if the current
        write operation fails.  Use at your own risk.
    ignore_divisions : bool, default False
        If False (default) raises error when previous divisions overlap with
        the new appended divisions. Ignored if append=False.
    partition_on : list, default None
        Construct directory-based partitioning by splitting on these fields'
        values. Each dask partition will result in one or more datafiles,
        there will be no global groupby.
    storage_options : dict, default None
        Key/value pairs to be passed on to the file-system backend, if any.
    custom_metadata : dict, default None
        Custom key/value metadata to include in all footer metadata (and
        in the global "_metadata" file, if applicable).  Note that the custom
        metadata may not contain the reserved b"pandas" key.
    write_metadata_file : bool or None, default None
        Whether to write the special ``_metadata`` file. If ``None`` (the
        default), a ``_metadata`` file will only be written if ``append=True``
        and the dataset already has a ``_metadata`` file.
    compute : bool, default True
        If ``True`` (default) then the result is computed immediately. If
        ``False`` then a ``dask.dataframe.Scalar`` object is returned for
        future computation.
    compute_kwargs : dict, default True
        Options to be passed in to the compute method
    schema : pyarrow.Schema, dict, "infer", or None, default "infer"
        Global schema to use for the output dataset. Defaults to "infer", which
        will infer the schema from the dask dataframe metadata. This is usually
        sufficient for common schemas, but notably will fail for ``object``
        dtype columns that contain things other than strings. These columns
        will require an explicit schema be specified. The schema for a subset
        of columns can be overridden by passing in a dict of column names to
        pyarrow types (for example ``schema={"field": pa.string()}``); columns
        not present in this dict will still be automatically inferred.
        Alternatively, a full ``pyarrow.Schema`` may be passed, in which case
        no schema inference will be done. Passing in ``schema=None`` will
        disable the use of a global file schema - each written file may use a
        different schema dependent on the dtypes of the corresponding
        partition.
    name_function : callable, default None
        Function to generate the filename for each output partition.
        The function should accept an integer (partition index) as input and
        return a string which will be used as the filename for the corresponding
        partition. Should preserve the lexicographic order of partitions.
        If not specified, files will created using the convention
        ``part.0.parquet``, ``part.1.parquet``, ``part.2.parquet``, ...
        and so on for each partition in the DataFrame.
    filesystem: "fsspec", "arrow", or fsspec.AbstractFileSystem backend to use.
    **kwargs :
        Extra options to be passed on to the specific backend.

    Examples
    --------
    >>> df = dd.read_csv(...)  # doctest: +SKIP
    >>> df.to_parquet('/path/to/output/', ...)  # doctest: +SKIP

    By default, files will be created in the specified output directory using the
    convention ``part.0.parquet``, ``part.1.parquet``, ``part.2.parquet``, ... and so on for
    each partition in the DataFrame. To customize the names of each file, you can use the
    ``name_function=`` keyword argument. The function passed to ``name_function`` will be
    used to generate the filename for each partition and should expect a partition's index
    integer as input and return a string which will be used as the filename for the corresponding
    partition. Strings produced by ``name_function`` must preserve the order of their respective
    partition indices.

    For example:

    >>> name_function = lambda x: f"data-{x}.parquet"
    >>> df.to_parquet('/path/to/output/', name_function=name_function)  # doctest: +SKIP

    will result in the following files being created::

        /path/to/output/
            ├── data-0.parquet
            ├── data-1.parquet
            ├── data-2.parquet
            └── ...

    See Also
    --------
    read_parquet: Read parquet data to dask.dataframe
    """
    from dask.dataframe.dask_expr._collection import new_collection

    engine = _set_parquet_engine(engine=engine, meta=df._meta)
    compute_kwargs = compute_kwargs or {}

    partition_on = partition_on or []
    if isinstance(partition_on, str):
        partition_on = [partition_on]

    if set(partition_on) - set(df.columns):
        raise ValueError(
            "Partitioning on non-existent column. "
            "partition_on=%s ."
            "columns=%s" % (str(partition_on), str(list(df.columns)))
        )

    if df.columns.inferred_type not in {"string", "empty"}:
        raise ValueError("parquet doesn't support non-string column names")

    if isinstance(engine, str):
        engine = get_engine(engine)

    if hasattr(path, "name"):
        path = stringify_path(path)

    fs, _paths, _, _ = engine.extract_filesystem(
        path,
        filesystem=filesystem,
        dataset_options={},
        open_file_options={},
        storage_options=storage_options,
    )
    assert len(_paths) == 1, "only one path"
    path = _paths[0]

    if overwrite:
        if append:
            raise ValueError("Cannot use both `overwrite=True` and `append=True`!")

        if fs.exists(path) and fs.isdir(path):
            # Check for any previous parquet ops reading from a file in the
            # output directory, since deleting those files now would result in
            # errors or incorrect results.
            for read_op in df.expr.find_operations(ReadParquet):
                read_path_with_slash = str(read_op.path).rstrip("/") + "/"
                write_path_with_slash = path.rstrip("/") + "/"
                if read_path_with_slash.startswith(write_path_with_slash):
                    raise ValueError(
                        "Cannot overwrite a path that you are reading "
                        "from in the same task graph."
                    )

            # Don't remove the directory if it's the current working directory
            if _is_local_fs(fs):
                working_dir = fs.expand_path(".")[0]
                if path.rstrip("/") == working_dir.rstrip("/"):
                    raise ValueError(
                        "Cannot clear the contents of the current working directory!"
                    )

            # It's safe to clear the output directory
            fs.rm(path, recursive=True)

        # Clear read_parquet caches in case we are
        # also reading from the overwritten path
        _cached_plan.clear()

    # Always skip divisions checks if divisions are unknown
    if not df.known_divisions:
        ignore_divisions = True

    # Save divisions and corresponding index name. This is necessary,
    # because we may be resetting the index to write the file
    division_info = {"divisions": df.divisions, "name": df.index.name}
    if division_info["name"] is None:
        # As of 0.24.2, pandas will rename an index with name=None
        # when df.reset_index() is called.  The default name is "index",
        # but dask will always change the name to the NONE_LABEL constant
        if NONE_LABEL not in df.columns:
            division_info["name"] = NONE_LABEL
        elif write_index:
            raise ValueError(
                "Index must have a name if __null_dask_index__ is a column."
            )
        else:
            warnings.warn(
                "If read back by Dask, column named __null_dask_index__ "
                "will be set to the index (and renamed to None)."
            )

    # There are some "reserved" names that may be used as the default column
    # name after resetting the index. However, we don't want to treat it as
    # a "special" name if the string is already used as a "real" column name.
    reserved_names = []
    for name in ["index", "level_0"]:
        if name not in df.columns:
            reserved_names.append(name)

    # If write_index==True (default), reset the index and record the
    # name of the original index in `index_cols` (we will set the name
    # to the NONE_LABEL constant if it is originally `None`).
    # `pyarrow` will revert the `reset_index` call
    # below if `index_cols` is populated (because pyarrow will want to handle
    # index preservation itself).  The column index
    # will be written to "pandas metadata" if write_index=True
    index_cols = []
    if write_index:
        real_cols = set(df.columns)
        none_index = list(df._meta.index.names) == [None]
        df = df.reset_index()
        if none_index:
            rename_columns = {c: NONE_LABEL for c in df.columns if c in reserved_names}
            df = df.rename(columns=rename_columns)
        index_cols = [c for c in set(df.columns) - real_cols]
    else:
        # Not writing index - might as well drop it
        df = df.reset_index(drop=True)

    if custom_metadata and b"pandas" in custom_metadata.keys():
        raise ValueError(
            "User-defined key/value metadata (custom_metadata) can not "
            "contain a b'pandas' key.  This key is reserved by Pandas, "
            "and overwriting the corresponding value can render the "
            "entire dataset unreadable."
        )

    # Engine-specific initialization steps to write the dataset.
    # Possibly create parquet metadata, and load existing stuff if appending
    i_offset, fmd, metadata_file_exists, extra_write_kwargs = engine.initialize_write(
        df,
        fs,
        path,
        append=append,
        ignore_divisions=ignore_divisions,
        partition_on=partition_on,
        division_info=division_info,
        index_cols=index_cols,
        schema=schema,
        custom_metadata=custom_metadata,
        **kwargs,
    )

    # By default we only write a metadata file when appending if one already
    # exists
    if append and write_metadata_file is None:
        write_metadata_file = metadata_file_exists

    # Check that custom name_function is valid,
    # and that it will produce unique names
    if name_function is not None:
        if not callable(name_function):
            raise ValueError("``name_function`` must be a callable with one argument.")
        filenames = [name_function(i + i_offset) for i in range(df.npartitions)]
        if len(set(filenames)) < len(filenames):
            raise ValueError("``name_function`` must produce unique filenames.")

    # If we are using a remote filesystem and retries is not set, bump it
    # to be more fault tolerant, as transient transport errors can occur.
    # The specific number 5 isn't hugely motivated: it's less than ten and more
    # than two.
    annotations = dask.config.get("annotations", {})
    if "retries" not in annotations and not _is_local_fs(fs):
        ctx = dask.annotate(retries=5)
    else:
        ctx = contextlib.nullcontext()

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Dask annotations ", category=UserWarning
        )
        with ctx:
            out = new_collection(
                ToParquet(
                    df,
                    path,
                    fs,
                    fmd,
                    engine,
                    i_offset,
                    partition_on,
                    write_metadata_file,
                    name_function,
                    toolz.merge(
                        kwargs,
                        {
                            "compression": compression,
                            "custom_metadata": custom_metadata,
                        },
                        extra_write_kwargs,
                    ),
                    append,
                )
            )

    if compute:
        out = out.compute(**compute_kwargs)

    # Invalidate the filesystem listing cache for the output path after write.
    # We do this before returning, even if `compute=False`. This helps ensure
    # that reading files that were just written succeeds.
    fs.invalidate_cache(path)

    return out


def _determine_type_mapper(
    *, user_types_mapper, dtype_backend, pyarrow_strings_enabled
):
    type_mappers = []

    def pyarrow_type_mapper(pyarrow_dtype):
        # Special case pyarrow strings to use more feature complete dtype
        # See https://github.com/pandas-dev/pandas/issues/50074
        if pyarrow_dtype == pa.string():
            return pd.StringDtype("pyarrow")
        else:
            return pd.ArrowDtype(pyarrow_dtype)

    # always use the user-defined mapper first, if available
    if user_types_mapper is not None:
        type_mappers.append(user_types_mapper)

    # next in priority is converting strings
    if pyarrow_strings_enabled:
        type_mappers.append({pa.string(): pd.StringDtype("pyarrow")}.get)
        type_mappers.append({pa.date32(): pd.ArrowDtype(pa.date32())}.get)
        type_mappers.append({pa.date64(): pd.ArrowDtype(pa.date64())}.get)

        def _convert_decimal_type(type):
            if pa.types.is_decimal(type):
                return pd.ArrowDtype(type)
            return None

        type_mappers.append(_convert_decimal_type)

    # and then nullable types
    if dtype_backend == "numpy_nullable":
        type_mappers.append(PYARROW_NULLABLE_DTYPE_MAPPING.get)
    elif dtype_backend == "pyarrow":
        type_mappers.append(pyarrow_type_mapper)

    def default_types_mapper(pyarrow_dtype):
        """Try all type mappers in order, starting from the user type mapper."""
        for type_converter in type_mappers:
            converted_type = type_converter(pyarrow_dtype)
            if converted_type is not None:
                return converted_type

    if len(type_mappers) > 0:
        return default_types_mapper


class ReadParquet(PartitionsFiltered, BlockwiseIO):
    _pickle_functools_cache = False
    _absorb_projections = True
    _filter_passthrough = False

    def _filter_passthrough_available(self, parent, dependents):
        return (
            super()._filter_passthrough_available(parent, dependents)
            and (isinstance(parent.predicate, (LE, GE, LT, GT, EQ, NE, And, Or)))
            and _DNF.extract_pq_filters(self, parent.predicate)._filters is not None
        )

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Index):
            # Column projection
            columns = determine_column_projection(self, parent, dependents)
            columns = [col for col in self.columns if col in columns]
            if set(columns) == set(self.columns):
                return
            return Index(
                self.substitute_parameters({"columns": columns, "_series": False})
            )

        if isinstance(parent, Projection):
            return super()._simplify_up(parent, dependents)

        if isinstance(parent, Filter) and self._filter_passthrough_available(
            parent, dependents
        ):
            # Predicate pushdown
            filters = _DNF.extract_pq_filters(self, parent.predicate)
            if filters._filters is not None:
                return self.substitute_parameters(
                    {
                        "filters": filters.combine(
                            self.operand("filters")
                        ).to_list_tuple()
                    }
                )

        if isinstance(parent, Lengths):
            _lengths = self._get_lengths()
            if _lengths:
                return Literal(_lengths)

        if isinstance(parent, Len):
            _lengths = self._get_lengths()
            if _lengths:
                return Literal(sum(_lengths))

    @property
    def columns(self):
        columns_operand = self.operand("columns")
        if columns_operand is None:
            return list(self._meta.columns)
        else:
            return _convert_to_list(columns_operand)

    @cached_property
    def _funcname(self):
        return "read_parquet"

    def __dask_tokenize__(self):
        if not self._determ_token:
            # TODO: Is there an actual need to overwrite this?
            self._determ_token = _tokenize_deterministic(
                funcname(type(self)), self.checksum, *self.operands[:-1]
            )
        return self._determ_token

    @cached_property
    def _name(self):
        return self._funcname + "-" + self.deterministic_token

    @property
    def checksum(self):
        return self._dataset_info["checksum"]

    def _tree_repr_argument_construction(self, i, op, header):
        if self._parameters[i] == "_dataset_info_cache":
            # Don't print this, very ugly
            return header
        return super()._tree_repr_argument_construction(i, op, header)

    @cached_property
    def _meta(self):
        meta = self._dataset_info["base_meta"]
        columns = _convert_to_list(self.operand("columns"))
        if self._series:
            assert len(columns) > 0
            return meta[columns[0]]
        elif columns is not None:
            return meta[columns]
        return meta

    @abstractmethod
    def _divisions(self):
        raise NotImplementedError

    @property
    def _fusion_compression_factor(self):
        if self.operand("columns") is None:
            return 1
        nr_original_columns = max(len(self._dataset_info["schema"].names) - 1, 1)
        return max(
            len(_convert_to_list(self.operand("columns"))) / nr_original_columns, 0.001
        )


class ReadParquetPyarrowFS(ReadParquet):
    _parameters = [
        "path",
        "columns",
        "filters",
        "categories",
        "index",
        "storage_options",
        "filesystem",
        "ignore_metadata_file",
        "calculate_divisions",
        "arrow_to_pandas",
        "pyarrow_strings_enabled",
        "kwargs",
        "_partitions",
        "_series",
        "_dataset_info_cache",
    ]
    _defaults = {
        "columns": None,
        "filters": None,
        "categories": None,
        "index": None,
        "storage_options": None,
        "filesystem": None,
        "ignore_metadata_file": True,
        "calculate_divisions": False,
        "arrow_to_pandas": None,
        "pyarrow_strings_enabled": True,
        "kwargs": None,
        "_partitions": None,
        "_series": False,
        "_dataset_info_cache": None,
    }
    _absorb_projections = True
    _filter_passthrough = True

    @cached_property
    def normalized_path(self):
        return _normalize_and_strip_protocol(self.path)

    @cached_property
    def fs(self):
        fs_input = self.operand("filesystem")
        if isinstance(fs_input, pa.fs.FileSystem):
            return fs_input
        else:
            fs = pa_fs.FileSystem.from_uri(self.path)[0]
            if storage_options := self.storage_options:
                # Use inferred region as the default
                region = {} if "region" in storage_options else {"region": fs.region}
                fs = type(fs)(**region, **storage_options)
            return fs

    def approx_statistics(self) -> dict:
        """Return an approximation of a single files statistics.

        This is determined by sampling a few files and averaging their statistics.

        Fields
        ------
        num_rows: avg
        num_row_groups: avg
        serialized_size: avg
        columns: list
            A list of all column statistics where individual fields are also
            averaged.


        Example
        -------
        {
            'num_rows': 1991129,
            'num_row_groups': 2.3333333333333335,
            'serialized_size': 6256.666666666667,
            'total_byte_size': 118030095,
            'columns': [
                {'total_compressed_size': 6284162.333333333,
                'total_uncompressed_size': 6347380.333333333,
                'path_in_schema': 'l_orderkey'},
                {'total_compressed_size': 9423516.333333334,
                'total_uncompressed_size': 9423063.333333334,
                'path_in_schema': 'l_partkey'},
                {'total_compressed_size': 9405796.666666666,
                'total_uncompressed_size': 9405346.666666666,
                'path_in_schema': 'l_suppkey'},
                ...
            ]
        }

        Returns
        -------
        dict
        """
        idxs = self.sample_statistics()
        files_to_consider = np.array(self._dataset_info["all_files"])[idxs]
        stats = [_STATS_CACHE[tokenize(finfo)] for finfo in files_to_consider]
        return _combine_stats(stats)

    def load_statistics(self, files=None, fragments=None):
        if files is None:
            files = self._dataset_info["all_files"]
        if fragments is None:
            fragments = self.fragments_unsorted
        # Collecting code samples is actually a little expensive (~100ms) and
        # we'd like this thing to be as low overhead as possible
        with dask.config.set({"distributed.diagnostics.computations.nframes": 0}):
            token_stats = flatten(
                dask.compute(_collect_statistics_plan(files, fragments))
            )
        for token, stats in token_stats:
            _STATS_CACHE[token] = stats

    def sample_statistics(self, n=3):
        """Sample statistics from the dataset.

        Sample N file statistics from the dataset. The files are chosen by
        sorting all files based on their binary file size and picking
        equidistant sampling points.

        In the special case of n=3 this corresponds to min/median/max.

        Returns
        -------
        ixs: list[int]
            The indices of files that were sampled
        """
        frags = self.fragments_unsorted
        finfos = np.array(self._dataset_info["all_files"])
        getsize = np.frompyfunc(lambda x: x.size, nin=1, nout=1)
        finfo_size_arr = getsize(finfos)
        finfo_argsort = finfo_size_arr.argsort()
        nfrags = len(frags)
        stepsize = max(nfrags // n, 1)
        finfos_sampled = []
        frags_samples = []
        ixs = []
        for i in range(0, nfrags, stepsize):
            sort_ix = finfo_argsort[i]
            # TODO: This is crude but the most conservative estimate
            sort_ix = sort_ix if sort_ix < nfrags else 0
            ixs.append(sort_ix)
            finfos_sampled.append(finfos[sort_ix])
            frags_samples.append(frags[sort_ix])
        self.load_statistics(finfos_sampled, frags_samples)
        return ixs

    @cached_property
    def raw_statistics(self):
        """Parquet statstics for every file in the dataset.
        The statistics do not include all the metadata that is stored in the
        file but only a subset. See also `_extract_stats`.
        """
        self.load_statistics()
        return [
            _STATS_CACHE[tokenize(finfo)] for finfo in self._dataset_info["all_files"]
        ]

    @cached_property
    def aggregated_statistics(self):
        """Aggregate statistics for every partition in the dataset.

        These statistics aggregated the row group statistics to partition level
        such that min/max/total_compressed_size/etc. corresponds to the entire
        partition instead of individual row groups.
        """
        return _aggregate_statistics_to_file(self.raw_statistics)

    def _get_lengths(self):
        # TODO: Filters that only filter partition_expr can be used as well
        if not self.filters:
            return tuple(stats["num_rows"] for stats in self.aggregated_statistics)

    @cached_property
    def _dataset_info(self):
        if rv := self.operand("_dataset_info_cache"):
            return rv
        dataset_info = {}

        path_normalized = self.normalized_path
        # We'll first treat the path as if it was a directory since this is the
        # most common case. Only if this fails, we'll treat it as a file. This
        # way, the happy path performs one remote request instead of two if we
        # were to check the type of the path first.
        try:
            # At this point we will post a listbucket request which includes the
            # same data as a HEAD request. The information included here (see
            # pyarrow FileInfo) are size, type, path and modified since
            # timestamps This isn't free but relatively cheap (200-300ms or less
            # for ~1k files)
            all_files = []
            for path in path_normalized:
                dataset_selector = pa_fs.FileSelector(path, recursive=True)
                all_files.extend(
                    [
                        finfo
                        for finfo in self.fs.get_file_info(dataset_selector)
                        if finfo.type == pa.fs.FileType.File
                    ]
                )
        except (NotADirectoryError, FileNotFoundError):
            all_files = [self.fs.get_file_info(path) for path in path_normalized]
        # TODO: At this point we could verify if we're dealing with a very
        # inhomogeneous datasets already without reading any further data

        metadata_file = False
        checksum = None
        dataset = None
        if not self.ignore_metadata_file:
            all_files = sorted(
                all_files, key=lambda x: x.base_name.endswith("_metadata")
            )
            if all_files[-1].base_name.endswith("_metadata"):
                metadata_file = all_files.pop()
                checksum = tokenize(metadata_file)
                dataset = pa_ds.parquet_dataset(
                    metadata_file.path,
                    filesystem=self.fs,
                )
                dataset_info["using_metadata_file"] = True
                dataset_info["fragments"] = [
                    FragmentWrapper(frag, self.fs) for frag in dataset.get_fragments()
                ]
                dataset_info["file_sizes"] = [None] * len(dataset_info["fragments"])

        if checksum is None:
            checksum = tokenize(all_files)
            dataset_info["file_sizes"] = [fi.size for fi in all_files]
        dataset_info["checksum"] = checksum
        if dataset is None:
            import pyarrow.parquet as pq

            dataset = pq.ParquetDataset(
                # TODO Just pass all_files once
                # https://github.com/apache/arrow/pull/40143 is available to
                # reduce latency
                [fi.path for fi in all_files],
                filesystem=self.fs,
                filters=self.filters,
            )
            dataset_info["using_metadata_file"] = False
            dataset_info["fragments"] = [
                FragmentWrapper(frag, self.fs) for frag in dataset.fragments
            ]
            dataset_info["all_files"] = all_files

        dataset_info["schema"] = dataset.schema
        dataset_info["base_meta"] = dataset.schema.empty_table().to_pandas()
        self._dataset_info_cache = dataset_info
        return dataset_info

    @cached_property
    def _division_from_stats(self):
        """If enabled, compute the divisions from the collected statistics.
        If divisions are possible to set, the second argument will be the
        argsort of the fragments such that the divisions are correct.

        Returns
        -------
        divisions
        argsort
        """
        if self.calculate_divisions and self.index is not None:
            index_name = self.index.name
            return _divisions_from_statistics(self.aggregated_statistics, index_name)
        return tuple([None] * (len(self.fragments_unsorted) + 1)), None

    def all_statistics_known(self) -> bool:
        """Whether all statistics have been fetched from remote store"""
        return all(
            tokenize(finfo) in _STATS_CACHE for finfo in self._dataset_info["all_files"]
        )

    def _fragment_sort_index(self):
        return self._division_from_stats[1]

    def _divisions(self):
        return self._division_from_stats[0]

    def _tune_up(self, parent):
        if self._fusion_compression_factor >= 1:
            return
        if isinstance(parent, FusedParquetIO):
            return
        return parent.substitute(self, FusedParquetIO(self))

    @cached_property
    def fragments(self):
        """Return all fragments in the dataset after filtering in the order as
        expected by the divisions.

        See also
        --------
        ReadParquetPyarrowFS.fragments_unsorted
        """
        if self._fragment_sort_index() is not None:
            return self.fragments_unsorted[self._fragment_sort_index()]
        return self.fragments_unsorted

    @property
    def fragments_unsorted(self):
        """All fragments in the dataset after filtering.

        No guarantees on ordering. This is ordered as the files are listed.

        See also
        --------
        ReadParquetPyarrowFS.fragments
        """
        if self.filters is not None:
            filter_expression = pq.filters_to_expression(self.filters)
            frags = [frag.fragment for frag in self._dataset_info["fragments"]]
            ds = pyarrow.dataset.FileSystemDataset(
                frags,
                self._dataset_info["schema"],
                format=frags[0].format,
                filesystem=self.fs,
            )
            return np.array(list(ds.get_fragments(filter=filter_expression)))
        return np.array([frag.fragment for frag in self._dataset_info["fragments"]])

    @property
    def _fusion_compression_factor(self):
        approx_stats = self.approx_statistics()
        total_uncompressed = 0
        after_projection = 0
        col_op = self.operand("columns") or self.columns
        for col in approx_stats["columns"]:
            total_uncompressed += col["total_uncompressed_size"]
            if col["path_in_schema"] in col_op:
                after_projection += col["total_uncompressed_size"]

        min_size = parse_bytes(
            dask.config.get("dataframe.parquet.minimum-partition-size")
        )
        total_uncompressed = max(total_uncompressed, min_size)
        return max(after_projection / total_uncompressed, 0.001)

    def _filtered_task(self, name: Key, index: int) -> Task:
        columns = self.columns.copy()
        index_name = self.index.name
        if self.index is not None:
            index_name = self.index.name
        schema = self._dataset_info["schema"].remove_metadata()
        if index_name:
            if columns is None:
                columns = list(schema.names)
            columns.append(index_name)
        # Note: We use kwargs here to have an easier time to take this apart in
        # FusedParquetIO
        return Task(
            name,
            ReadParquetPyarrowFS._table_to_pandas,
            Task(
                None,
                ReadParquetPyarrowFS._fragment_to_table,
                fragment_wrapper=FragmentWrapper(
                    self.fragments[index], filesystem=self.fs
                ),
                filters=self.filters,
                columns=columns,
                schema=schema,
            ),
            index_name=index_name,
            arrow_to_pandas=self.arrow_to_pandas,
            dtype_backend=self.kwargs.get("dtype_backend"),
            pyarrow_strings_enabled=self.pyarrow_strings_enabled,
            _data_producer=True,
        )

    @staticmethod
    def _fragment_to_table(fragment_wrapper, filters, columns, schema):
        _maybe_adjust_cpu_count()
        if isinstance(fragment_wrapper, FragmentWrapper):
            fragment = fragment_wrapper.fragment
        else:
            fragment = fragment_wrapper
        if isinstance(filters, list):
            filters = pq.filters_to_expression(filters)
        return fragment.to_table(
            schema=schema,
            columns=columns,
            filter=filters,
            # Batch size determines how many rows are read at once and will
            # cause the underlying array to be split into chunks of this size
            # (max). We'd like to avoid fragmentation as much as possible and
            # and to set this to something like inf but we have to set a finite,
            # positive number.
            # In the presence of row groups, the underlying array will still be
            # chunked per rowgroup
            batch_size=10_000_000,
            fragment_scan_options=pa.dataset.ParquetFragmentScanOptions(
                pre_buffer=True,
                cache_options=pa.CacheOptions(
                    hole_size_limit=parse_bytes("4 MiB"),
                    range_size_limit=parse_bytes("32.00 MiB"),
                ),
            ),
            # TODO: Reconsider this. The OMP_NUM_THREAD variable makes it harmful to enable this
            use_threads=True,
        )

    @staticmethod
    def _table_to_pandas(
        table, index_name, arrow_to_pandas, dtype_backend, pyarrow_strings_enabled
    ):
        if arrow_to_pandas is None:
            arrow_to_pandas = {}
        else:
            arrow_to_pandas = arrow_to_pandas.copy()
        # This can mess up index setting, etc.
        arrow_to_pandas.pop("ignore_metadata", None)
        df = table.to_pandas(
            types_mapper=_determine_type_mapper(
                user_types_mapper=arrow_to_pandas.pop("types_mapper", None),
                dtype_backend=dtype_backend,
                pyarrow_strings_enabled=pyarrow_strings_enabled,
            ),
            use_threads=arrow_to_pandas.get("use_threads", False),
            self_destruct=arrow_to_pandas.get("self_destruct", True),
            **arrow_to_pandas,
            ignore_metadata=True,
        )
        if index_name is not None:
            df = df.set_index(index_name)
        return df


class ReadParquetFSSpec(ReadParquet):
    """Read a parquet dataset"""

    _parameters = [
        "path",
        "columns",
        "filters",
        "categories",
        "index",
        "storage_options",
        "calculate_divisions",
        "ignore_metadata_file",
        "metadata_task_size",
        "split_row_groups",
        "blocksize",
        "aggregate_files",
        "parquet_file_extension",
        "filesystem",
        "engine",
        "kwargs",
        "_partitions",
        "_series",
        "_dataset_info_cache",
        "_pq_length_stats",
    ]
    _defaults = {
        "columns": None,
        "filters": None,
        "categories": None,
        "index": None,
        "storage_options": None,
        "calculate_divisions": False,
        "ignore_metadata_file": False,
        "metadata_task_size": None,
        "split_row_groups": "infer",
        "blocksize": "default",
        "aggregate_files": None,
        "parquet_file_extension": (".parq", ".parquet", ".pq"),
        "filesystem": "fsspec",
        "engine": "pyarrow",
        "kwargs": {"dtype_backend": None},
        "_partitions": None,
        "_series": False,
        "_dataset_info_cache": None,
        "_pq_length_stats": None,
    }

    @property
    def engine(self):
        _engine = self.operand("engine")
        if isinstance(_engine, str):
            return get_engine(_engine)
        return _engine

    def _divisions(self):
        return self._plan["divisions"]

    @property
    def _dataset_info(self):
        if rv := self.operand("_dataset_info_cache"):
            return rv
        # Process and split user options
        (
            dataset_options,
            read_options,
            open_file_options,
            other_options,
        ) = _split_user_options(**(self.kwargs or {}))

        # Extract global filesystem and paths
        fs, paths, dataset_options, open_file_options = self.engine.extract_filesystem(
            self.path,
            self.filesystem,
            dataset_options,
            open_file_options,
            self.storage_options,
        )
        read_options["open_file_options"] = open_file_options
        paths = sorted(paths, key=natural_sort_key)  # numeric rather than glob ordering

        auto_index_allowed = False
        index_operand = self.operand("index")
        if index_operand is None:
            # User is allowing auto-detected index
            auto_index_allowed = True
        if index_operand and isinstance(index_operand, str):
            index = [index_operand]
        else:
            index = index_operand

        blocksize = self.blocksize
        if self.split_row_groups in ("infer", "adaptive"):
            # Using blocksize to plan partitioning
            if self.blocksize == "default":
                if hasattr(self.engine, "default_blocksize"):
                    blocksize = self.engine.default_blocksize()
                else:
                    blocksize = "128MiB"
        else:
            # Not using blocksize - Set to `None`
            blocksize = None

        # Collect general dataset info
        args = (
            paths,
            fs,
            self.categories,
            index,
            self.calculate_divisions,
            self.filters,
            self.split_row_groups,
            blocksize,
            self.aggregate_files,
            self.ignore_metadata_file,
            self.metadata_task_size,
            self.parquet_file_extension,
            {
                "read": read_options,
                "dataset": dataset_options,
                **other_options,
            },
        )
        dataset_info = self.engine._collect_dataset_info(*args)
        checksum = []
        files_for_checksum = []
        if dataset_info["has_metadata_file"]:
            if isinstance(self.path, list):
                files_for_checksum = [
                    next(path for path in self.path if path.endswith("_metadata"))
                ]
            else:
                files_for_checksum = [self.path + fs.sep + "_metadata"]
        else:
            files_for_checksum = dataset_info["ds"].files

        for file in files_for_checksum:
            # The checksum / file info is usually already cached by the fsspec
            # FileSystem dir_cache since this info was already asked for in
            # _collect_dataset_info
            checksum.append(fs.checksum(file))
        dataset_info["checksum"] = tokenize(checksum)

        # Infer meta, accounting for index and columns arguments.
        meta = self.engine._create_dd_meta(dataset_info)
        index = dataset_info["index"]
        index = [index] if isinstance(index, str) else index
        meta, index, all_columns = set_index_columns(
            meta, index, None, auto_index_allowed
        )
        if meta.index.name == NONE_LABEL:
            meta.index.name = None
        dataset_info["base_meta"] = meta
        dataset_info["index"] = index
        dataset_info["all_columns"] = all_columns
        dataset_info["calculate_divisions"] = self.calculate_divisions

        dataset_token = tokenize(dataset_info)
        if dataset_token not in _cached_plan:
            parts, stats, common_kwargs = self.engine._construct_collection_plan(
                dataset_info
            )

            # Make sure parts and stats are aligned
            parts, stats = _align_statistics(parts, stats)

            # Use statistics to aggregate partitions
            parts, stats = _aggregate_row_groups(parts, stats, dataset_info)

            # Drop filtered partitions (aligns with `dask.dataframe` behavior)
            if self.filters and stats:
                parts, stats = apply_filters(parts, stats, self.filters)

            # Use statistics to calculate divisions
            divisions = _calculate_divisions(stats, dataset_info, len(parts))

            empty = False
            if len(divisions) < 2:
                # empty dataframe - just use meta
                divisions = (None, None)
                parts = [meta]
                empty = True

            _control_cached_plan(dataset_token)
            _cached_plan[dataset_token] = {
                "empty": empty,
                "parts": parts,
                "statistics": stats,
                "divisions": divisions,
                "common_kwargs": common_kwargs,
            }
        dataset_info["plan"] = _cached_plan[dataset_token]
        self._dataset_info_cache = dataset_info
        return dataset_info

    def _filtered_task(self, name: Key, index: int) -> Task:
        tsk = Task(name, self._io_func, self._plan["parts"][index], _data_producer=True)
        if self._series:
            return Task(name, operator.getitem, tsk, self.columns[0])
        return tsk

    @property
    def _io_func(self):
        if self._plan["empty"]:
            return identity
        dataset_info = self._dataset_info
        return ParquetFunctionWrapper(
            self.engine,
            dataset_info["fs"],
            dataset_info["base_meta"],
            self.columns,
            dataset_info["index"],
            dataset_info["kwargs"]["dtype_backend"],
            {},  # All kwargs should now be in `common_kwargs`
            self._plan["common_kwargs"],
        )

    @property
    def _plan(self):
        return self._dataset_info["plan"]

    def _get_lengths(self) -> tuple | None:
        """Return known partition lengths using parquet statistics"""
        if not self.filters:
            return tuple(
                length
                for i, length in enumerate(self._pq_length_stats)
                if not self._filtered or i in self._partitions
            )
        return None

    @cached_property
    def _pq_length_stats(self):
        """Ensure that partition-length statistics are up to date"""

        if self._plan["statistics"]:
            # Already have statistics from original API call
            return tuple(
                stat["num-rows"]
                for i, stat in enumerate(self._plan["statistics"])
                if not self._filtered or i in self._partitions
            )
        else:
            # Need to go back and collect statistics
            return tuple(stat["num-rows"] for stat in _collect_pq_statistics(self))


#
# Helper functions
#


def _set_parquet_engine(engine=None, meta=None):
    # Use `engine` or `meta` input to set the parquet engine
    if engine == "fastparquet":
        raise NotImplementedError("Fastparquet engine is not supported")

    if engine is None:
        if (
            meta is not None and typename(meta).split(".")[0] == "cudf"
        ) or dask.config.get("dataframe.backend", "pandas") == "cudf":
            from dask_cudf.io.parquet import CudfEngine

            engine = CudfEngine
        else:
            engine = "pyarrow"
    return engine


def _align_statistics(parts, statistics):
    # Make sure parts and statistics are aligned
    # (if statistics is not empty)
    if statistics and len(parts) != len(statistics):
        statistics = []
    if statistics:
        result = list(
            zip(
                *[
                    (part, stats)
                    for part, stats in zip(parts, statistics)
                    if stats["num-rows"] > 0
                ]
            )
        )
        parts, statistics = result or [[], []]
    return parts, statistics


def _aggregate_row_groups(parts, statistics, dataset_info):
    # Aggregate parts/statistics if we are splitting by row-group
    blocksize = (
        dataset_info["blocksize"] if dataset_info["split_row_groups"] is True else None
    )
    split_row_groups = dataset_info["split_row_groups"]
    fs = dataset_info["fs"]
    aggregation_depth = dataset_info["aggregation_depth"]

    if statistics:
        if blocksize or (split_row_groups and int(split_row_groups) > 1):
            parts, statistics = aggregate_row_groups(
                parts, statistics, blocksize, split_row_groups, fs, aggregation_depth
            )
    return parts, statistics


def _calculate_divisions(statistics, dataset_info, npartitions):
    # Use statistics to define divisions
    divisions = None
    if statistics and dataset_info.get("gather_statistics", False):
        calculate_divisions = dataset_info.get("calculate_divisions", None)
        index = dataset_info["index"]
        process_columns = index if index and len(index) == 1 else None
        if (calculate_divisions is not False) and process_columns:
            for sorted_column_info in sorted_columns(
                statistics, columns=process_columns
            ):
                if sorted_column_info["name"] in index:
                    divisions = sorted_column_info["divisions"]
                    break

    return divisions or (None,) * (npartitions + 1)


#
# Filtering logic
#


class _DNF:
    """Manage filters in Disjunctive Normal Form (DNF)"""

    class _Or(frozenset):
        """Frozen set of disjunctions"""

        def to_list_tuple(self) -> list:
            # DNF "or" is List[List[Tuple]]
            def _maybe_list(val):
                if isinstance(val, tuple) and val and isinstance(val[0], (tuple, list)):
                    return list(val)
                return [val]

            return [
                (
                    _maybe_list(val.to_list_tuple())
                    if hasattr(val, "to_list_tuple")
                    else _maybe_list(val)
                )
                for val in self
            ]

    class _And(frozenset):
        """Frozen set of conjunctions"""

        def to_list_tuple(self) -> list:
            # DNF "and" is List[Tuple]
            return tuple(  # type: ignore
                val.to_list_tuple() if hasattr(val, "to_list_tuple") else val
                for val in self
            )

    _filters: _And | _Or | None  # Underlying filter expression

    def __init__(self, filters: _And | _Or | list | tuple | None) -> None:
        self._filters = self.normalize(filters)

    def to_list_tuple(self) -> list:
        return self._filters.to_list_tuple()  # type: ignore

    def __bool__(self) -> bool:
        return bool(self._filters)

    @classmethod
    def normalize(cls, filters: _And | _Or | list | tuple | None):
        """Convert raw filters to the `_Or(_And)` DNF representation"""
        if not filters:
            result = None
        elif isinstance(filters, list):
            conjunctions = filters if isinstance(filters[0], list) else [filters]
            result = cls._Or([cls._And(conjunction) for conjunction in conjunctions])
        elif isinstance(filters, tuple):
            if isinstance(filters[0], tuple):
                raise TypeError("filters must be List[Tuple] or List[List[Tuple]]")
            result = cls._Or((cls._And((filters,)),))
        elif isinstance(filters, cls._Or):
            result = cls._Or(se for e in filters for se in cls.normalize(e))
        elif isinstance(filters, cls._And):
            total = []
            for c in itertools.product(*[cls.normalize(e) for e in filters]):
                total.append(cls._And(se for e in c for se in e))
            result = cls._Or(total)
        else:
            raise TypeError(f"{type(filters)} not a supported type for _DNF")
        return result

    def combine(self, other: _DNF | _And | _Or | list | tuple | None) -> _DNF:
        """Combine with another _DNF object"""
        if not isinstance(other, _DNF):
            other = _DNF(other)
        assert isinstance(other, _DNF)
        if self._filters is None:
            result = other._filters
        elif other._filters is None:
            result = self._filters
        else:
            result = self._And([self._filters, other._filters])
        return _DNF(result)

    @classmethod
    def extract_pq_filters(cls, pq_expr: ReadParquet, predicate_expr: Expr) -> _DNF:
        _filters = None
        if isinstance(predicate_expr, (LE, GE, LT, GT, EQ, NE)):
            if (
                not isinstance(predicate_expr.right, Expr)
                and isinstance(predicate_expr.left, Projection)
                and predicate_expr.left.frame._name == pq_expr._name
            ):
                op = predicate_expr._operator_repr
                column = predicate_expr.left.columns[0]
                value = predicate_expr.right
                _filters = (column, op, value)
            elif (
                not isinstance(predicate_expr.left, Expr)
                and isinstance(predicate_expr.left, Projection)
                and predicate_expr.left.frame._name == pq_expr._name
            ):
                # Simple dict to make sure field comes first in filter
                flip = {LE: GE, LT: GT, GE: LE, GT: LT}
                op = predicate_expr  # type: ignore
                op = flip.get(op, op)._operator_repr  # type: ignore
                column = predicate_expr.right.columns[0]
                value = predicate_expr.left
                _filters = (column, op, value)

        elif isinstance(predicate_expr, (And, Or)):
            left = cls.extract_pq_filters(pq_expr, predicate_expr.left)._filters
            right = cls.extract_pq_filters(pq_expr, predicate_expr.right)._filters
            if left and right:
                if isinstance(predicate_expr, And):
                    _filters = cls._And([left, right])  # type: ignore
                else:
                    _filters = cls._Or([left, right])  # type: ignore

        return _DNF(_filters)


#
# Parquet-statistics handling
#


def _collect_pq_statistics(
    expr: ReadParquet, columns: list | None = None
) -> list[dict] | None:
    """Collect Parquet statistic for dataset paths"""

    # Be strict about columns argument
    if columns:
        if not isinstance(columns, list):
            raise ValueError(f"Expected columns to be a list, got {type(columns)}.")
        allowed = {expr._meta.index.name} | set(expr.columns)
        if not set(columns).issubset(allowed):
            raise ValueError(f"columns={columns} must be a subset of {allowed}")

    if expr._plan["empty"]:
        return []

    # Collect statistics using layer information
    fs = expr._io_func.fs
    parts = [
        part
        for i, part in enumerate(expr._plan["parts"])
        if not expr._filtered or i in expr._partitions
    ]

    # Execute with delayed for large and remote datasets
    parallel = int(False if _is_local_fs(fs) else 16)
    if parallel:
        # Group parts corresponding to the same file.
        # A single task should always parse statistics
        # for all these parts at once (since they will
        # all be in the same footer)
        groups = defaultdict(list)
        for part in parts:
            for p in [part] if isinstance(part, dict) else part:
                path = p.get("piece")[0]  # type: ignore[index]
                groups[path].append(p)
        group_keys = list(groups.keys())

        # Compute and return flattened result
        func = delayed(_read_partition_stats_group)
        result = dask.compute(
            [
                func(
                    list(
                        itertools.chain(
                            *[groups[k] for k in group_keys[i : i + parallel]]
                        )
                    ),
                    fs,
                    columns=columns,
                )
                for i in range(0, len(group_keys), parallel)
            ]
        )[0]
        return list(itertools.chain(*result))
    else:
        # Serial computation on client
        return _read_partition_stats_group(parts, fs, columns=columns)


def _read_partition_stats_group(parts, fs, columns=None):
    """Parse the statistics for a group of files"""

    def _read_partition_stats(part, fs, columns=None):
        # Helper function to read Parquet-metadata
        # statistics for a single partition

        if not isinstance(part, list):
            part = [part]

        column_stats = {}
        num_rows = 0
        columns = columns or []
        for p in part:
            piece = p["piece"]
            path = piece[0]
            row_groups = None if piece[1] == [None] else piece[1]
            with fs.open(path, default_cache="none") as f:
                md = pq.ParquetFile(f).metadata
            if row_groups is None:
                row_groups = list(range(md.num_row_groups))
            for rg in row_groups:
                row_group = md.row_group(rg)
                num_rows += row_group.num_rows
                for i in range(row_group.num_columns):
                    col = row_group.column(i)
                    name = col.path_in_schema
                    if name in columns:
                        if col.statistics and col.statistics.has_min_max:
                            if name in column_stats:
                                column_stats[name]["min"] = min(
                                    column_stats[name]["min"], col.statistics.min
                                )
                                column_stats[name]["max"] = max(
                                    column_stats[name]["max"], col.statistics.max
                                )
                            else:
                                column_stats[name] = {
                                    "min": col.statistics.min,
                                    "max": col.statistics.max,
                                }

        # Convert dict-of-dict to list-of-dict to be consistent
        # with current `dd.read_parquet` convention (for now)
        column_stats_list = [
            {
                "name": name,
                "min": column_stats[name]["min"],
                "max": column_stats[name]["max"],
            }
            for name in column_stats.keys()
        ]
        return {"num-rows": num_rows, "columns": column_stats_list}

    # Helper function used by _extract_statistics
    return [_read_partition_stats(part, fs, columns=columns) for part in parts]


def _normalize_and_strip_protocol(path):
    if not isinstance(path, (list, tuple)):
        path = [path]

    result = []
    for p in path:
        protocol_separators = ["://", "::"]
        for sep in protocol_separators:
            split = p.split(sep, 1)
            if len(split) > 1:
                p = split[1]
                break
        result.append(p.rstrip("/"))
    return result


def _divisions_from_statistics(aggregated_stats, index_name):
    col_ix = -1
    peak_rg = aggregated_stats[0]
    for ix, col in enumerate(peak_rg["columns"]):
        if col["path_in_schema"] == index_name:
            col_ix = ix
            break
    else:
        raise ValueError(
            f"Index column {index_name} not found in statistics"  # noqa: E713
        )
    last_max = None
    minmax = []
    for file_stats in aggregated_stats:
        file_min = file_stats["columns"][col_ix]["statistics"]["min"]
        file_max = file_stats["columns"][col_ix]["statistics"]["max"]

        minmax.append((file_min, file_max))
    divisions = []
    minmax = pd.Series(minmax)

    argsort = minmax.argsort()
    sorted_minmax = minmax[argsort]
    if not sorted_minmax.is_monotonic_increasing:
        return tuple([None] * (len(aggregated_stats) + 1)), None
    for file_min, file_max in sorted_minmax:
        divisions.append(file_min)
        last_max = file_max
    divisions.append(last_max)
    return tuple(divisions), argsort


def _extract_stats(original):
    """Take the raw file statistics as returned by pyarrow (as a dict) and
    filter it to what we care about. The full stats are a bit too verbose and we
    don't need all of it."""
    # TODO: dicts are pretty memory inefficient. Move to dataclass?
    file_level_stats = ["num_rows", "num_row_groups", "serialized_size"]
    rg_stats = [
        "num_rows",
        "total_byte_size",
        "sorting_columns",
    ]
    col_meta = [
        "num_values",
        "total_compressed_size",
        "total_uncompressed_size",
        "path_in_schema",
    ]
    col_stats = [
        "min",
        "max",
        "null_count",
        "num_values",
        "distinct_count",
    ]

    out = {}
    for name in file_level_stats:
        out[name] = original[name]
    out["row_groups"] = rgs = []
    for rg in original["row_groups"]:
        rg_out = {}
        rgs.append(rg_out)
        for name in rg_stats:
            rg_out[name] = rg[name]
        rg_out["columns"] = []
        for col in rg["columns"]:
            col_out = {}
            rg_out["columns"].append(col_out)
            for name in col_meta:
                col_out[name] = col[name]
            col_out["statistics"] = {}
            if col["statistics"] is None:
                continue
            for name in col_stats:
                col_out["statistics"][name] = col["statistics"][name]

    return out


def _agg_dicts(dicts, agg_funcs):
    result = {}
    for d in dicts:
        for k, v in d.items():
            if k not in result:
                result[k] = [v]
            else:
                result[k].append(v)
    result2 = {}
    for k, v in result.items():
        agg = agg_funcs.get(k)
        if agg:
            result2[k] = agg(v)
    return result2


def _aggregate_columns(cols, agg_cols):
    combine = []
    i = 0
    while True:
        inner = []
        combine.append(inner)
        try:
            for col in cols:
                inner.append(col[i])
        except IndexError:
            combine.pop()
            break
        i += 1
    return [_agg_dicts(c, agg_cols) for c in combine]


def _get_min_max_value(x, func):
    x = [y for y in x if y is not None]
    return func(x) if len(x) > 0 else None


def _aggregate_statistics_to_file(stats):
    """Aggregate RG information to file level."""

    agg_stats = {
        "min": lambda x: _get_min_max_value(x, min),
        "max": lambda x: _get_min_max_value(x, max),
    }
    agg_cols = {
        "total_compressed_size": sum,
        "total_uncompressed_size": sum,
        "statistics": partial(_agg_dicts, agg_funcs=agg_stats),
        "path_in_schema": lambda x: set(x).pop(),
    }
    agg_func = {
        "num_rows": sum,
        "total_byte_size": sum,
        "columns": partial(_aggregate_columns, agg_cols=agg_cols),
    }
    aggregated_stats = []
    for file_stat in stats:
        file_stat = file_stat.copy()
        aggregated_stats.append(file_stat)

        file_stat.update(_agg_dicts(file_stat.pop("row_groups"), agg_func))
    return aggregated_stats


@dask.delayed
def _gather_statistics(frags):
    @dask.delayed
    def _collect_statistics(token_fragment):
        if isinstance(token_fragment[1], FragmentWrapper):
            token_fragment[1] = token_fragment[1].fragment
        return token_fragment[0], _extract_stats(token_fragment[1].metadata.to_dict())

    return dask.compute(
        list(_collect_statistics(frag) for frag in frags), scheduler="threading"
    )[0]


def _collect_statistics_plan(file_infos, fragments):
    """Collect statistics for a list of files and their corresponding fragments"""

    to_collect = []
    for finfo, frag in zip(file_infos, fragments):
        if (token := tokenize(finfo)) not in _STATS_CACHE:
            to_collect.append((token, frag))
    return [
        _gather_statistics(batch)
        for batch in toolz.itertoolz.partition_all(20, to_collect)
    ]


def _combine_stats(stats):
    """Combine multiple file-level statistics into a single dict of metrics that
    represent the average values of the parquet statistics"""
    agg_cols = {
        "total_compressed_size": statistics.mean,
        "total_uncompressed_size": statistics.mean,
        "path_in_schema": lambda x: set(x).pop(),
    }
    return _agg_dicts(
        _aggregate_statistics_to_file(stats),
        {
            "num_rows": statistics.mean,
            "num_row_groups": statistics.mean,
            "serialized_size": statistics.mean,
            "total_byte_size": statistics.mean,
            "columns": partial(_aggregate_columns, agg_cols=agg_cols),
        },
    )
