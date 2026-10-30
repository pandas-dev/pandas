"""parquet compat"""

from __future__ import annotations

from collections import abc
import io
import itertools
import json
import os
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Self,
    overload,
)
import warnings
from warnings import (
    catch_warnings,
    filterwarnings,
)

from pandas._config.config import _global_config as config

from pandas._libs import lib
from pandas.compat._optional import import_optional_dependency
from pandas.errors import (
    AbstractMethodError,
    Pandas4Warning,
)
from pandas.util._decorators import set_module
from pandas.util._exceptions import find_stack_level
from pandas.util._validators import check_dtype_backend

from pandas import (
    DataFrame,
    RangeIndex,
)

from pandas.io._util import (
    arrow_table_to_pandas,
    suppress_pyarrow_values_warning,
)
from pandas.io.common import (
    IOHandles,
    check_parent_directory,
    get_handle,
    is_fsspec_url,
    is_url,
    stringify_path,
)
from pandas.io.parsers.readers import validate_integer

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Generator,
        Hashable,
        Iterable,
    )

    import pyarrow

    from pandas._typing import (
        DtypeBackend,
        FilePath,
        ParquetCompressionOptions,
        ReadBuffer,
        StorageOptions,
        WriteBuffer,
    )


def get_engine(engine: str) -> BaseImpl:
    """return our implementation"""
    if engine == "auto":
        engine = config["io"]["parquet"]["engine"]

    if engine == "auto":
        # try engines in this order
        engine_classes = [PyArrowImpl, FastParquetImpl]

        error_msgs = ""
        for engine_class in engine_classes:
            try:
                return engine_class()
            except ImportError as err:
                error_msgs += "\n - " + str(err)

        raise ImportError(
            "Unable to find a usable engine; "
            "tried using: 'pyarrow', 'fastparquet'.\n"
            "A suitable version of "
            "pyarrow or fastparquet is required for parquet "
            "support.\n"
            "Trying to import the above resulted in these errors:"
            f"{error_msgs}"
        )

    if engine == "pyarrow":
        return PyArrowImpl()
    elif engine == "fastparquet":
        warnings.warn(
            "The 'fastparquet' engine is deprecated and will be removed in a "
            "future version. Use 'pyarrow' instead.",
            Pandas4Warning,
            stacklevel=find_stack_level(),
        )
        return FastParquetImpl()

    raise ValueError("engine must be one of 'pyarrow', 'fastparquet'")


def _get_path_or_handle(
    path: FilePath | ReadBuffer[bytes] | WriteBuffer[bytes],
    fs: Any,
    storage_options: StorageOptions | None = None,
    mode: str = "rb",
    is_dir: bool = False,
) -> tuple[
    FilePath | ReadBuffer[bytes] | WriteBuffer[bytes], IOHandles[bytes] | None, Any
]:
    """File handling for PyArrow."""
    path_or_handle = stringify_path(path)
    if fs is not None:
        pa_fs = import_optional_dependency("pyarrow.fs", errors="ignore")
        fsspec = import_optional_dependency("fsspec", errors="ignore")
        if pa_fs is not None and isinstance(fs, pa_fs.FileSystem):
            if storage_options:
                raise NotImplementedError(
                    "storage_options not supported with a pyarrow FileSystem."
                )
        elif fsspec is not None and isinstance(fs, fsspec.spec.AbstractFileSystem):
            pass
        else:
            raise ValueError(
                f"filesystem must be a pyarrow or fsspec FileSystem, "
                f"not a {type(fs).__name__}"
            )
    if is_fsspec_url(path_or_handle) and fs is None:
        pa_error = None
        if storage_options is None:
            pa = import_optional_dependency("pyarrow")
            pa_fs = import_optional_dependency("pyarrow.fs")

            try:
                fs, path_or_handle = pa_fs.FileSystem.from_uri(path)
            except (TypeError, pa.ArrowException):
                pass
            except OSError as err:
                # Only "hdfs:///path" (no host) resolves differently in fsspec
                # (GH#58078); otherwise this is a real error, e.g. no libhdfs
                if not str(path_or_handle).startswith("hdfs:///"):
                    raise
                pa_error = err
        if fs is None:
            try:
                fsspec = import_optional_dependency("fsspec")
                fs, path_or_handle = fsspec.core.url_to_fs(
                    path_or_handle, **(storage_options or {})
                )
            except Exception as err:
                if pa_error is None:
                    raise
                # keep pyarrow's error, e.g. a missing JVM
                raise err from pa_error
    elif storage_options and (not is_url(path_or_handle) or mode != "rb"):
        # can't write to a remote url
        # without making use of fsspec at the moment
        raise ValueError("storage_options passed with buffer, or non-supported URL")

    handles = None
    if (
        not fs
        and not is_dir
        and isinstance(path_or_handle, str)
        and not os.path.isdir(path_or_handle)
    ):
        if is_url(path_or_handle):
            # pyarrow cannot read non-fsspec URLs (e.g. http/https), so let
            # get_handle download them into a buffer for pyarrow to consume.
            handles = get_handle(
                path_or_handle, mode, is_text=False, storage_options=storage_options
            )
            path_or_handle = handles.handle
        else:
            # Local path: keep the I/O in pyarrow's C++ layer rather than the
            # Python one (GH#47702), and open the file only once -- adding a
            # get_handle open on top would clobber pyarrow's data to 0 bytes
            # (GH#65810). get_handle would also expand "~" and check the parent
            # directory on write, so reproduce both below.
            path_or_handle = os.path.expanduser(path_or_handle)
            if "w" in mode or "a" in mode or "x" in mode:
                check_parent_directory(path_or_handle)
                # Open the destination instead of handing over its path:
                # write_table deletes a path-like target when the write raises,
                # destroying a pre-existing file it never managed to open
                # (GH#69022). pa.OSFile opens it in C++ rather than through
                # builtins.open, so this stays a single native open.
                pa = import_optional_dependency("pyarrow")
                stream = pa.OSFile(path_or_handle, mode)
                handles = IOHandles(
                    handle=stream,
                    compression={"method": None},
                    created_handles=[stream],
                )
                path_or_handle = stream
    return path_or_handle, handles, fs


def _is_dir(path: str, fs: Any) -> bool:
    """Whether ``path`` is a directory on ``fs``, or locally if ``fs`` is None."""
    if fs is None:
        return os.path.isdir(path)
    pa_fs = import_optional_dependency("pyarrow.fs")
    if isinstance(fs, pa_fs.FileSystem):
        return fs.get_file_info(path).type == pa_fs.FileType.Directory
    # otherwise an fsspec filesystem, as checked in _get_path_or_handle
    return fs.isdir(path)


def _fragment_batches(
    fragments: Iterable[Any], batch_size: int, **scanner_kwargs: Any
) -> Generator[pyarrow.Table]:
    """
    Tables of ``batch_size`` rows, the last possibly shorter, read from each
    row group of each dataset fragment in turn.

    A fragment's scan reads row groups ahead of its consumer, with or without
    threads and however small ``batch_readahead`` is, so a slow consumer would
    end up holding most of the file. Scanning each row group on its own holds
    about one at a time, like ``ParquetFile.iter_batches()``. A scan's batches
    end with its row group, so they are joined across row groups and
    fragments into tables of ``batch_size`` rows, as ``iter_batches()`` joins
    them across the row groups of a file.

    Each row group is read through a ``RecordBatchReader`` that is closed when
    it is exhausted or when this generator is closed, which releases its file
    right away rather than once the scan is garbage collected.
    """
    pa = import_optional_dependency("pyarrow")
    pending: list[pyarrow.RecordBatch] = []
    pending_rows = 0
    for fragment in fragments:
        for row_group in fragment.split_by_row_group():
            scanner = row_group.scanner(batch_size=batch_size, **scanner_kwargs)
            with scanner.to_reader() as reader:
                for batch in reader:
                    while pending_rows + batch.num_rows >= batch_size:
                        split = batch_size - pending_rows
                        yield pa.Table.from_batches([*pending, batch.slice(0, split)])
                        batch = batch.slice(split)
                        pending, pending_rows = [], 0
                    # an empty slice would still hold on to its row group
                    if batch.num_rows:
                        pending.append(batch)
                        pending_rows += batch.num_rows
    if pending_rows:
        yield pa.Table.from_batches(pending)


def _batches_or_empty(
    batches: Generator[pyarrow.RecordBatch | pyarrow.Table],
    read_empty: Callable[[], pyarrow.Table],
) -> Generator[pyarrow.RecordBatch | pyarrow.Table]:
    """
    ``batches``, or if there are none the empty table ``read_empty()`` reads.

    Like ``read_csv(chunksize=...)``, empty data then gives one empty chunk
    with the columns of a non-chunked read rather than no chunks at all.
    Closing this generator closes ``batches``.
    """
    empty = True
    try:
        for batch in batches:
            empty = False
            yield batch
    finally:
        batches.close()
    if empty:
        yield read_empty()


def _full_read_range_index(
    schema: pyarrow.Schema,
    num_rows: Callable[[], int],
    to_pandas_kwargs: dict[str, Any] | None,
) -> tuple[int, int, Hashable]:
    """
    Start, step and name of the RangeIndex a non-chunked read would restore.

    ``to_pandas`` applies a RangeIndex stored in the pandas metadata only when
    it is as long as the table being converted, so a single batch gets a
    default 0-based RangeIndex instead.
    """
    metadata = schema.pandas_metadata
    if metadata and not (to_pandas_kwargs or {}).get("ignore_metadata"):
        index_columns = metadata["index_columns"]
        if len(index_columns) == 1 and isinstance(index_columns[0], dict):
            descr = index_columns[0]
            stored = range(descr["start"], descr["stop"], descr["step"])
            if descr["kind"] == "range" and len(stored) == num_rows():
                return descr["start"], descr["step"], descr["name"]
    return 0, 1, None


def _arrow_to_frame(
    table: pyarrow.Table | pyarrow.RecordBatch,
    dtype_backend: DtypeBackend | lib.NoDefault,
    to_pandas_kwargs: dict[str, Any] | None,
) -> DataFrame:
    """Convert a pyarrow table or batch to a DataFrame, restoring ``attrs``."""
    df_metadata = None
    if table.schema.metadata:
        if b"PANDAS_ATTRS" in table.schema.metadata:
            df_metadata = table.schema.metadata[b"PANDAS_ATTRS"]

    with catch_warnings():
        filterwarnings(
            "ignore",
            "make_block is deprecated",
            Pandas4Warning,
        )
        result = arrow_table_to_pandas(
            table,
            dtype_backend=dtype_backend,
            to_pandas_kwargs=to_pandas_kwargs,
        )

    if df_metadata is not None:
        result.attrs = json.loads(df_metadata)
    return result


class BaseImpl:
    @staticmethod
    def validate_dataframe(df: DataFrame) -> None:
        if not isinstance(df, DataFrame):
            raise ValueError("to_parquet only supports IO with DataFrames")

    def write(self, df: DataFrame, path, compression, **kwargs) -> None:
        raise AbstractMethodError(self)

    def read(self, path, columns=None, **kwargs) -> DataFrame:
        raise AbstractMethodError(self)


@set_module("pandas.api.typing")
class ParquetFileReader(abc.Iterator[DataFrame]):
    """
    Iterator over DataFrame chunks read from a parquet file or dataset.

    A ``ParquetFileReader`` is returned by :func:`~pandas.read_parquet` when
    ``chunksize`` is passed; it is not usually instantiated directly. It
    wraps pyarrow's ``ParquetFile.iter_batches()`` (single file) or the
    fragments of a ``ParquetDataset`` (partitioned directory) and
    converts each batch to a DataFrame lazily, one at a time, closing the
    underlying file(s) once the iterator is exhausted or used as a context
    manager.

    Parameters
    ----------
    batch_iter : Iterator[pyarrow.RecordBatch | pyarrow.Table]
        Iterator over pyarrow record batches or tables, or an empty table for
        empty data, from either
        ``pyarrow.parquet.ParquetFile.iter_batches()`` or the fragments of a
        ``pyarrow.parquet.ParquetDataset``.
    source : pyarrow.parquet.ParquetFile or generator
        The object ``batch_iter`` reads from, closed along with ``handles``.
    dtype_backend : {'numpy_nullable', 'pyarrow'} or lib.no_default
        Back-end data type applied to each yielded DataFrame. See
        :func:`~pandas.read_parquet` for details.
    to_pandas_kwargs : dict or None
        Keyword arguments forwarded to ``pyarrow.Table.to_pandas`` when
        converting each batch.
    handles : IOHandles or None
        Open file handle(s) to close once the iterator is exhausted or
        used as a context manager.
    index_range : tuple of (int, int, Hashable)
        Start, step and name of the RangeIndex a non-chunked read would
        restore; chunks with a RangeIndex continue it.

    See Also
    --------
    read_parquet : Load a parquet object from the file path.

    Notes
    -----
    Only returned when ``engine="pyarrow"``; :func:`~pandas.read_parquet`
    raises ``NotImplementedError`` for ``chunksize`` with
    ``engine="fastparquet"``, together with ``filters``, or for a single file
    with ``columns`` selecting a field inside a nested column or with the
    ``schema``, ``partitioning`` or ``ignore_prefixes`` keywords.

    Examples
    --------
    >>> from io import BytesIO
    >>> df = pd.DataFrame({"a": range(10)})
    >>> buf = BytesIO(df.to_parquet())
    >>> reader = pd.read_parquet(buf, chunksize=5)
    >>> for chunk in reader:
    ...     print(len(chunk))
    5
    5
    """

    def __init__(
        self,
        batch_iter,
        source,
        dtype_backend: DtypeBackend | lib.NoDefault,
        to_pandas_kwargs: dict[str, Any] | None,
        handles: IOHandles[bytes] | None,
        index_range: tuple[int, int, Hashable],
    ) -> None:
        self._batch_iter = batch_iter
        self._source = source
        self._dtype_backend: DtypeBackend | lib.NoDefault = dtype_backend
        self._to_pandas_kwargs = to_pandas_kwargs
        self._handles = handles
        self._index_range = index_range
        self._closed = False
        self._rows_read = 0

    def __next__(self) -> DataFrame:
        try:
            batch = next(self._batch_iter)
        except StopIteration:
            self.close()
            raise
        df = _arrow_to_frame(batch, self._dtype_backend, self._to_pandas_kwargs)
        # A batch gets a RangeIndex starting at 0 rather than its part of the
        # stored one (see _full_read_range_index), so continue the full read's
        # RangeIndex across chunks, like read_csv(chunksize=...) does.
        if isinstance(df.index, RangeIndex):
            start, step, name = self._index_range
            first = start + self._rows_read * step
            df.index = RangeIndex(first, first + len(df) * step, step, name=name)
        self._rows_read += len(df)
        return df

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *args) -> None:
        self.close()

    def close(self) -> None:
        """
        Close the underlying file(s), if any are open.

        Called automatically when the iterator is exhausted or when the
        ``ParquetFileReader`` is used as a context manager.

        See Also
        --------
        read_parquet : Load a parquet object from the file path.

        Examples
        --------
        >>> reader = pd.read_parquet("data.parquet", chunksize=5)  # doctest: +SKIP
        >>> reader.close()  # doctest: +SKIP
        """
        if not self._closed:
            # pyarrow opens local and filesystem paths itself, in which case
            # there are no pandas handles and only this closes the file
            self._source.close()
            if self._handles is not None:
                self._handles.close()
        self._closed = True


class PyArrowImpl(BaseImpl):
    def __init__(self) -> None:
        import_optional_dependency(
            "pyarrow", extra="pyarrow is required for parquet support."
        )
        import pyarrow.parquet

        # import utils to register the pyarrow extension types
        import pandas.core.arrays.arrow.extension_types  # pyright: ignore[reportUnusedImport] # noqa: F401

        self.api = pyarrow

    def write(
        self,
        df: DataFrame,
        path: FilePath | WriteBuffer[bytes],
        compression: ParquetCompressionOptions = "snappy",
        index: bool | None = None,
        storage_options: StorageOptions | None = None,
        partition_cols: list[str] | None = None,
        filesystem=None,
        **kwargs,
    ) -> None:
        self.validate_dataframe(df)

        from_pandas_kwargs: dict[str, Any] = {"schema": kwargs.pop("schema", None)}
        if index is not None:
            from_pandas_kwargs["preserve_index"] = index

        with suppress_pyarrow_values_warning():
            table = self.api.Table.from_pandas(df, **from_pandas_kwargs)

        if df.attrs:
            df_metadata = {"PANDAS_ATTRS": json.dumps(df.attrs)}
            existing_metadata = table.schema.metadata
            merged_metadata = {**existing_metadata, **df_metadata}
            table = table.replace_schema_metadata(merged_metadata)

        if partition_cols is None and kwargs:
            # pyarrow opens the destination before validating kwargs, so a
            # misspelled kwarg would clobber an existing file or be masked by
            # an error about the path (GH#45815). Validate against a buffer,
            # skipping kwargs that are single-use (encryption_properties on
            # pyarrow<20) or record each write (metadata_collector).
            self.api.parquet.write_table(
                table.schema.empty_table(),
                self.api.BufferOutputStream(),
                compression=compression,
                **{**kwargs, "metadata_collector": None, "encryption_properties": None},
            )

        path_or_handle, handles, filesystem = _get_path_or_handle(
            path,
            filesystem,
            storage_options=storage_options,
            mode="wb",
            is_dir=partition_cols is not None,
        )
        try:
            if partition_cols is not None:
                # writes to multiple files under the given path
                self.api.parquet.write_to_dataset(
                    table,
                    path_or_handle,
                    compression=compression,
                    partition_cols=partition_cols,
                    filesystem=filesystem,
                    **kwargs,
                )
            else:
                # write to single output file
                self.api.parquet.write_table(
                    table,
                    path_or_handle,
                    compression=compression,
                    filesystem=filesystem,
                    **kwargs,
                )
        finally:
            if handles is not None:
                handles.close()

    @overload
    def read(
        self,
        path,
        columns=None,
        filters=None,
        dtype_backend: DtypeBackend | lib.NoDefault = ...,
        storage_options: StorageOptions | None = ...,
        filesystem=None,
        to_pandas_kwargs: dict[str, Any] | None = ...,
        *,
        chunksize: None = ...,
        **kwargs,
    ) -> DataFrame: ...

    @overload
    def read(
        self,
        path,
        columns=None,
        filters=None,
        dtype_backend: DtypeBackend | lib.NoDefault = ...,
        storage_options: StorageOptions | None = ...,
        filesystem=None,
        to_pandas_kwargs: dict[str, Any] | None = ...,
        *,
        chunksize: int,
        **kwargs,
    ) -> ParquetFileReader: ...

    def read(
        self,
        path,
        columns=None,
        filters=None,
        dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
        storage_options: StorageOptions | None = None,
        filesystem=None,
        to_pandas_kwargs: dict[str, Any] | None = None,
        *,
        chunksize: int | None = None,
        **kwargs,
    ) -> DataFrame | ParquetFileReader:
        if chunksize is not None and filters is not None:
            raise NotImplementedError("filters is not supported when chunksize is set")
        kwargs["use_pandas_metadata"] = True

        path_or_handle, handles, filesystem = _get_path_or_handle(
            path,
            filesystem,
            storage_options=storage_options,
            mode="rb",
        )

        if chunksize is not None:
            # Of read_table()'s keywords only use_threads is about iterating;
            # the others configure reading, so they go where the data is opened.
            use_threads = kwargs.pop("use_threads", True)
            del kwargs["use_pandas_metadata"]
            # pyarrow's default pre-buffering keeps the row groups it buffers
            # until the file is closed, so the reader would end up holding the
            # whole file
            kwargs.setdefault("pre_buffer", False)
            try:
                if isinstance(path_or_handle, str) and _is_dir(
                    path_or_handle, filesystem
                ):
                    # unlike pyarrow.dataset.dataset(), ParquetDataset takes the
                    # keywords and defaults of read_table(), e.g. hive partitions
                    dataset = self.api.parquet.ParquetDataset(
                        path_or_handle, filesystem=filesystem, **kwargs
                    )
                    schema = dataset.schema
                    # to_batches() has no use_pandas_metadata, so add the index
                    # columns that read_table() reads along with columns=
                    pandas_metadata = schema.pandas_metadata
                    if columns is not None and pandas_metadata:
                        columns = list(columns) + [
                            col
                            for col in pandas_metadata["index_columns"]
                            # a RangeIndex is stored as a dict, not as a column
                            if isinstance(col, str) and col not in columns
                        ]
                    fragments = dataset.fragments
                    index_range = _full_read_range_index(
                        schema,
                        lambda: sum(fragment.count_rows() for fragment in fragments),
                        to_pandas_kwargs,
                    )
                    source = batch_iter = _batches_or_empty(
                        _fragment_batches(
                            fragments,
                            schema=schema,
                            columns=columns,
                            batch_size=chunksize,
                            use_threads=use_threads,
                        ),
                        lambda: dataset.read(columns=columns, use_pandas_metadata=True),
                    )
                else:
                    # read_table() reads a file through a ParquetDataset, which
                    # takes these, but ParquetFile doesn't
                    dataset_only = [
                        key
                        for key in ("schema", "partitioning", "ignore_prefixes")
                        if key in kwargs
                    ]
                    if dataset_only:
                        raise NotImplementedError(
                            "keyword arguments are not supported when chunksize "
                            f"is set for a file: {dataset_only}"
                        )
                    parquet_file = self.api.parquet.ParquetFile(
                        path_or_handle, filesystem=filesystem, **kwargs
                    )
                    if columns is not None:
                        # iter_batches() leaves out names it doesn't know where
                        # read_table() raises, and gives the dotted path to a
                        # field inside a nested column as that column cut down
                        # to the field, where read_table() gives the field
                        known = {
                            prefix
                            for column in parquet_file.schema
                            for prefix in itertools.accumulate(
                                column.path.split("."), lambda a, b: f"{a}.{b}"
                            )
                        }
                        top_level = parquet_file.schema_arrow.names
                        missing = [col for col in columns if col not in known]
                        nested = [
                            col
                            for col in columns
                            if col in known and col not in top_level
                        ]
                        if missing or nested:
                            parquet_file.close()
                            if missing:
                                raise ValueError(
                                    f"Columns not found in the parquet file: {missing}"
                                )
                            raise NotImplementedError(
                                "selecting fields inside nested columns is not "
                                f"supported when chunksize is set for a file: {nested}"
                            )
                    batch_iter = _batches_or_empty(
                        parquet_file.iter_batches(
                            batch_size=chunksize,
                            columns=columns,
                            use_threads=use_threads,
                            use_pandas_metadata=True,
                        ),
                        lambda: parquet_file.read(
                            columns=columns, use_pandas_metadata=True
                        ),
                    )
                    index_range = _full_read_range_index(
                        parquet_file.schema_arrow,
                        lambda: parquet_file.metadata.num_rows,
                        to_pandas_kwargs,
                    )
                    source = parquet_file
            except Exception:
                if handles is not None:
                    handles.close()
                raise

            return ParquetFileReader(
                batch_iter,
                source=source,
                dtype_backend=dtype_backend,
                to_pandas_kwargs=to_pandas_kwargs,
                handles=handles,
                index_range=index_range,
            )

        try:
            pa_table = self.api.parquet.read_table(
                path_or_handle,
                columns=columns,
                filesystem=filesystem,
                filters=filters,
                **kwargs,
            )
            return _arrow_to_frame(pa_table, dtype_backend, to_pandas_kwargs)
        finally:
            if handles is not None:
                handles.close()


class FastParquetImpl(BaseImpl):
    def __init__(self) -> None:
        # since pandas is a dependency of fastparquet
        # we need to import on first use
        fastparquet = import_optional_dependency(
            "fastparquet", extra="fastparquet is required for parquet support."
        )
        self.api = fastparquet

    def write(
        self,
        df: DataFrame,
        path,
        compression: Literal["snappy", "gzip", "brotli"] | None = "snappy",
        index=None,
        partition_cols=None,
        storage_options: StorageOptions | None = None,
        filesystem=None,
        **kwargs,
    ) -> None:
        self.validate_dataframe(df)

        if "partition_on" in kwargs and partition_cols is not None:
            raise ValueError(
                "Cannot use both partition_on and "
                "partition_cols. Use partition_cols for partitioning data"
            )
        if "partition_on" in kwargs:
            partition_cols = kwargs.pop("partition_on")

        if partition_cols is not None:
            kwargs["file_scheme"] = "hive"

        if filesystem is not None:
            raise NotImplementedError(
                "filesystem is not implemented for the fastparquet engine."
            )

        # cannot use get_handle as write() does not accept file buffers
        path = stringify_path(path)
        if is_fsspec_url(path):
            fsspec = import_optional_dependency("fsspec")

            # if filesystem is provided by fsspec, file must be opened in 'wb' mode.
            kwargs["open_with"] = lambda path, _: fsspec.open(
                path, "wb", **(storage_options or {})
            ).open()
        elif storage_options:
            raise ValueError(
                "storage_options passed with file object or non-fsspec file path"
            )

        with catch_warnings(record=True):
            self.api.write(
                path,
                df,
                compression=compression,
                write_index=index,
                partition_on=partition_cols,
                **kwargs,
            )

    def read(
        self,
        path,
        columns=None,
        filters=None,
        storage_options: StorageOptions | None = None,
        filesystem=None,
        to_pandas_kwargs: dict[str, Any] | None = None,
        *,
        chunksize: int | None = None,
        **kwargs,
    ) -> DataFrame:
        if chunksize is not None:
            raise NotImplementedError(
                "chunksize is only supported with engine='pyarrow'"
            )
        parquet_kwargs: dict[str, Any] = {}
        dtype_backend = kwargs.pop("dtype_backend", lib.no_default)
        # We are disabling nullable dtypes for fastparquet pending discussion
        parquet_kwargs["pandas_nulls"] = False
        if dtype_backend is not lib.no_default:
            raise ValueError(
                "The 'dtype_backend' argument is not supported for the "
                "fastparquet engine"
            )
        if filesystem is not None:
            raise NotImplementedError(
                "filesystem is not implemented for the fastparquet engine."
            )
        if to_pandas_kwargs is not None:
            raise NotImplementedError(
                "to_pandas_kwargs is not implemented for the fastparquet engine."
            )
        path = stringify_path(path)
        handles = None
        if is_fsspec_url(path):
            fsspec = import_optional_dependency("fsspec")

            parquet_kwargs["fs"] = fsspec.open(path, "rb", **(storage_options or {})).fs
        elif isinstance(path, str) and not os.path.isdir(path):
            # use get_handle only when we are very certain that it is not a directory
            # fsspec resources can also point to directories
            # this branch is used for example when reading from non-fsspec URLs
            handles = get_handle(
                path, "rb", is_text=False, storage_options=storage_options
            )
            path = handles.handle

        try:
            parquet_file = self.api.ParquetFile(path, **parquet_kwargs)
            with catch_warnings():
                filterwarnings(
                    "ignore",
                    "make_block is deprecated",
                    Pandas4Warning,
                )
                return parquet_file.to_pandas(
                    columns=columns, filters=filters, **kwargs
                )
        finally:
            if handles is not None:
                handles.close()


def to_parquet(
    df: DataFrame,
    path: FilePath | WriteBuffer[bytes] | None = None,
    engine: str | lib.NoDefault = lib.no_default,
    compression: ParquetCompressionOptions = "snappy",
    index: bool | None = None,
    storage_options: StorageOptions | None = None,
    partition_cols: list[str] | None = None,
    filesystem: Any = None,
    **kwargs,
) -> bytes | None:
    """
    Write a DataFrame to the parquet format.

    Parameters
    ----------
    df : DataFrame
    path : str, path object, file-like object, or None, default None
        String, path object (implementing ``os.PathLike[str]``), or file-like
        object implementing a binary ``write()`` function. If None, the result
        is returned as bytes. If a string, it will be used as the root directory
        path when writing a partitioned dataset. The engine fastparquet does
        not accept file-like objects.

        The string could be a URL. Valid URL schemes include http, ftp, s3,
        gs, and file. For file URLs, a host is expected. A local file could be:
        ``file://localhost/path/to/table.parquet``. A remote example could be:
        ``s3://bucket/path/to/table.parquet``.

        Certain URL schemes may require additional packages. For example, S3
        URLs require the ``s3fs`` library. See
        :ref:`install.optional_dependencies` for a full list.
    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet library to use. If 'auto', then the option
        ``io.parquet.engine`` is used. The default ``io.parquet.engine``
        behavior is to try 'pyarrow', falling back to 'fastparquet' if
        'pyarrow' is unavailable.

        When using the ``'pyarrow'`` engine and no storage options are provided
        and a filesystem is implemented by both ``pyarrow.fs`` and ``fsspec``
        (e.g. "s3://"), then the ``pyarrow.fs`` filesystem is attempted first.
        Use the filesystem keyword with an instantiated fsspec filesystem
        if you wish to use its implementation.

        .. deprecated:: 3.1.0
            The ``'fastparquet'`` and ``'auto'`` engine options are
            deprecated. Use ``'pyarrow'`` or do not pass ``engine``
            to use the default.

    compression : {'snappy', 'gzip', 'brotli', 'lz4', 'zstd', None},
        default 'snappy'. Name of the compression to use. Use ``None``
        for no compression.
    index : bool, default None
        If ``True``, include the dataframe's index(es) in the file output. If
        ``False``, they will not be written to the file.
        If ``None``, similar to ``True`` the dataframe's index(es)
        will be saved. However, instead of being saved as values,
        the RangeIndex will be stored as a range in the metadata so it
        doesn't require much space and is faster. Other indexes will
        be included as columns in the file output.
    partition_cols : str or list, optional, default None
        Column names by which to partition the dataset.
        Columns are partitioned in the order they are given.
        Must be None if path is not a string.
    storage_options : dict, optional
        Extra options that make sense for a particular storage connection, e.g.
        host, port, username, password, etc. For HTTP(S) URLs the key-value
        pairs are forwarded to ``urllib.request.Request`` as header options.
        For other URLs (e.g. starting with "s3://", and "gcs://") the
        key-value pairs are forwarded to ``fsspec.open``. Please see ``fsspec``
        and ``urllib`` for more details, and for more examples on storage
        options refer `here <https://pandas.pydata.org/docs/user_guide/io.html?
        highlight=storage_options#reading-writing-remote-files>`_.
    filesystem : fsspec or pyarrow filesystem, default None
        Filesystem object to use when reading the parquet file. Only implemented
        for ``engine="pyarrow"``.

        .. versionadded:: 2.1.0

    **kwargs
        Additional keyword arguments passed to the engine:

        * For ``engine="pyarrow"``: passed to :func:`pyarrow.parquet.write_table`
          or :func:`pyarrow.parquet.write_to_dataset` (when using partition_cols)
        * For ``engine="fastparquet"``: passed to :func:`fastparquet.write`

    Returns
    -------
    bytes if no path argument is provided else None
    """
    if isinstance(partition_cols, str):
        partition_cols = [partition_cols]
    if engine is not lib.no_default and engine == "auto":
        warnings.warn(
            "engine='auto' is deprecated and will be removed in a "
            "future version. Pass engine='pyarrow' or do not pass "
            "'engine' to use the default.",
            Pandas4Warning,
            stacklevel=find_stack_level(),
        )
    if engine is lib.no_default:
        engine = "auto"
    impl = get_engine(engine)

    path_or_buf: FilePath | WriteBuffer[bytes] = io.BytesIO() if path is None else path

    impl.write(
        df,
        path_or_buf,
        compression=compression,
        index=index,
        partition_cols=partition_cols,
        storage_options=storage_options,
        filesystem=filesystem,
        **kwargs,
    )

    if path is None:
        assert isinstance(path_or_buf, io.BytesIO)
        return path_or_buf.getvalue()
    else:
        return None


@overload
def read_parquet(
    path: FilePath | ReadBuffer[bytes],
    engine: str | lib.NoDefault = ...,
    columns: list[str] | None = ...,
    storage_options: StorageOptions | None = ...,
    dtype_backend: DtypeBackend | lib.NoDefault = ...,
    filesystem: Any = ...,
    filters: list[tuple[Any, ...]] | list[list[tuple[Any, ...]]] | None = ...,
    to_pandas_kwargs: dict[str, Any] | None = ...,
    *,
    chunksize: None = ...,
    **kwargs,
) -> DataFrame: ...


@overload
def read_parquet(
    path: FilePath | ReadBuffer[bytes],
    engine: str | lib.NoDefault = ...,
    columns: list[str] | None = ...,
    storage_options: StorageOptions | None = ...,
    dtype_backend: DtypeBackend | lib.NoDefault = ...,
    filesystem: Any = ...,
    filters: list[tuple[Any, ...]] | list[list[tuple[Any, ...]]] | None = ...,
    to_pandas_kwargs: dict[str, Any] | None = ...,
    *,
    chunksize: int,
    **kwargs,
) -> ParquetFileReader: ...


@set_module("pandas")
def read_parquet(
    path: FilePath | ReadBuffer[bytes],
    engine: str | lib.NoDefault = lib.no_default,
    columns: list[str] | None = None,
    storage_options: StorageOptions | None = None,
    dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
    filesystem: Any = None,
    filters: list[tuple[Any, ...]] | list[list[tuple[Any, ...]]] | None = None,
    to_pandas_kwargs: dict[str, Any] | None = None,
    *,
    chunksize: int | None = None,
    **kwargs,
) -> DataFrame | ParquetFileReader:
    """
    Load a parquet object from the file path, returning a DataFrame.

    This function requires the `pyarrow <https://arrow.apache.org/docs/python/>`_
    library.

    The function automatically handles reading the data from a parquet file
    and creates a DataFrame with the appropriate structure.

    Parameters
    ----------
    path : str, path object or file-like object
        String, path object (implementing ``os.PathLike[str]``), or file-like
        object implementing a binary ``read()`` function.
        The string could be a URL. Valid URL schemes include http, ftp, s3,
        gs, and file. For file URLs, a host is expected. A local file could be:
        ``file://localhost/path/to/table.parquet``.

        Certain URL schemes may require additional packages. For example, S3
        URLs require the ``s3fs`` library. See
        :ref:`install.optional_dependencies` for a full list.

        A file URL can also be a path to a directory that contains multiple
        partitioned parquet files. Both pyarrow and fastparquet support
        paths to directories as well as file URLs. A directory path could be:
        ``file://localhost/path/to/tables`` or ``s3://bucket/partition_dir``.
    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet library to use. If 'auto', then the option
        ``io.parquet.engine`` is used. The default ``io.parquet.engine``
        behavior is to try 'pyarrow', falling back to 'fastparquet' if
        'pyarrow' is unavailable.

        When using the ``'pyarrow'`` engine and no storage options are provided
        and a filesystem is implemented by both ``pyarrow.fs`` and ``fsspec``
        (e.g. "s3://"), then the ``pyarrow.fs`` filesystem is attempted first.
        Use the filesystem keyword with an instantiated fsspec filesystem
        if you wish to use its implementation.

        .. deprecated:: 3.1.0
            The ``'fastparquet'`` and ``'auto'`` engine options are
            deprecated. Use ``'pyarrow'`` or do not pass ``engine``
            to use the default.

    columns : list, default=None
        If not None, only these columns will be read from the file.
    storage_options : dict, optional
        Extra options that make sense for a particular storage connection, e.g.
        host, port, username, password, etc. For HTTP(S) URLs the key-value
        pairs are forwarded to ``urllib.request.Request`` as header options.
        For other URLs (e.g. starting with "s3://", and "gcs://") the
        key-value pairs are forwarded to ``fsspec.open``. Please see ``fsspec``
        and ``urllib`` for more details, and for more examples on storage
        options refer `here <https://pandas.pydata.org/docs/user_guide/io.html?
        highlight=storage_options#reading-writing-remote-files>`_.
    dtype_backend : {'numpy_nullable', 'pyarrow'}
        Back-end data type applied to the resultant :class:`DataFrame`
        (still experimental). If not specified, the default behavior
        is to not use nullable data types. If specified, the behavior
        is as follows:

        * ``"numpy_nullable"``: returns nullable-dtype-backed :class:`DataFrame`
        * ``"pyarrow"``: returns pyarrow-backed nullable
          :class:`ArrowDtype` :class:`DataFrame`

        .. versionadded:: 2.0

    filesystem : fsspec or pyarrow filesystem, default None
        Filesystem object to use when reading the parquet file. Only implemented
        for ``engine="pyarrow"``.

        .. versionadded:: 2.1.0

    filters : List[Tuple] or List[List[Tuple]], default None
        To filter out data.
        Filter syntax: [[(column, op, val), ...],...]
        where op is [==, =, >, >=, <, <=, !=, in, not in]
        The innermost tuples are transposed into a set of filters applied
        through an `AND` operation.
        The outer list combines these sets of filters through an `OR`
        operation.
        A single list of tuples can also be used, meaning that no `OR`
        operation between set of filters is to be conducted.

        Using this argument will NOT result in row-wise filtering of the final
        partitions unless ``engine="pyarrow"`` is also specified.  For
        other engines, filtering is only performed at the partition level, that is,
        to prevent the loading of some row-groups and/or files.

        .. versionadded:: 2.1.0

    to_pandas_kwargs : dict | None, default None
        Keyword arguments to pass through to :func:`pyarrow.Table.to_pandas`
        when ``engine="pyarrow"``.

        .. versionadded:: 3.0.0

    chunksize : int, optional
        If specified, return a ``ParquetFileReader`` object that lazily
        reads and yields DataFrames of at most ``chunksize`` rows each,
        instead of loading the whole file into memory at once. Only
        supported with ``engine="pyarrow"``, and not supported together
        with ``filters``. Files are read without pre-buffering unless
        ``pre_buffer=True`` is passed, which can be faster on high-latency
        filesystems such as S3 but keeps the data read so far in memory.

        .. versionadded:: 3.1.0

    **kwargs
        Additional keyword arguments passed to the engine:

        * For ``engine="pyarrow"``: passed to :func:`pyarrow.parquet.read_table`.
          With ``chunksize``, ``use_threads`` is used for iterating and the
          others are passed to :class:`pyarrow.parquet.ParquetFile` for a file
          or :class:`pyarrow.parquet.ParquetDataset` for a directory;
          ``schema``, ``partitioning`` and ``ignore_prefixes`` are not
          supported for a file.
        * For ``engine="fastparquet"``: passed to
          :meth:`fastparquet.ParquetFile.to_pandas`

    Returns
    -------
    DataFrame or ParquetFileReader
        DataFrame based on parquet file, or an iterator of DataFrames if
        ``chunksize`` is specified.

    See Also
    --------
    DataFrame.to_parquet : Create a parquet object that serializes a DataFrame.

    Examples
    --------
    >>> original_df = pd.DataFrame({"foo": range(5), "bar": range(5, 10)})
    >>> original_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9
    >>> df_parquet_bytes = original_df.to_parquet()
    >>> from io import BytesIO
    >>> restored_df = pd.read_parquet(BytesIO(df_parquet_bytes))
    >>> restored_df
       foo  bar
    0    0    5
    1    1    6
    2    2    7
    3    3    8
    4    4    9
    >>> restored_df.equals(original_df)
    True
    >>> restored_bar = pd.read_parquet(BytesIO(df_parquet_bytes), columns=["bar"])
    >>> restored_bar
        bar
    0    5
    1    6
    2    7
    3    8
    4    9
    >>> restored_bar.equals(original_df[["bar"]])
    True

    The function uses `kwargs` that are passed directly to the engine.
    In the following example, we use the `filters` argument of the pyarrow
    engine to filter the rows of the DataFrame.

    Since `pyarrow` is the default engine, we can omit the `engine` argument.
    Note that the `filters` argument is implemented by the `pyarrow` engine,
    which can benefit from multithreading and also potentially be more
    economical in terms of memory.

    >>> sel = [("foo", ">", 2)]
    >>> restored_part = pd.read_parquet(BytesIO(df_parquet_bytes), filters=sel)
    >>> restored_part
        foo  bar
    0    3    8
    1    4    9
    """

    if engine is not lib.no_default and engine == "auto":
        warnings.warn(
            "engine='auto' is deprecated and will be removed in a "
            "future version. Pass engine='pyarrow' or do not pass "
            "'engine' to use the default.",
            Pandas4Warning,
            stacklevel=find_stack_level(),
        )
    if engine is lib.no_default:
        engine = "auto"
    impl = get_engine(engine)
    check_dtype_backend(dtype_backend)
    chunksize = validate_integer("chunksize", chunksize, 1)

    return impl.read(
        path,
        columns=columns,
        filters=filters,
        storage_options=storage_options,
        dtype_backend=dtype_backend,
        filesystem=filesystem,
        to_pandas_kwargs=to_pandas_kwargs,
        chunksize=chunksize,
        **kwargs,
    )
