""" parquet compat """

from typing import Any, AnyStr, Dict, List, Optional
from warnings import catch_warnings

from pandas._typing import FilePathOrBuffer
from pandas.compat._optional import import_optional_dependency
from pandas.errors import AbstractMethodError

from pandas import DataFrame, get_option

from pandas.io.common import _expand_user, get_filepath_or_buffer, is_fsspec_url


def get_engine(engine: str) -> "BaseImpl":
    """ return our implementation """
    if engine == "auto":
        engine = get_option("io.parquet.engine")

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
        return FastParquetImpl()

    raise ValueError("engine must be one of 'pyarrow', 'fastparquet'")


class BaseImpl:
    @staticmethod
    def validate_dataframe(df: DataFrame):

        if not isinstance(df, DataFrame):
            raise ValueError("to_parquet only supports IO with DataFrames")

        # must have value column names (strings only)
        if df.columns.inferred_type not in {"string", "empty"}:
            raise ValueError("parquet must have string column names")

        # index level names must be strings
        valid_names = all(
            isinstance(name, str) for name in df.index.names if name is not None
        )
        if not valid_names:
            raise ValueError("Index level names must be strings")

    def write(self, df: DataFrame, path, compression, **kwargs):
        raise AbstractMethodError(self)

    def read(self, path, columns=None, **kwargs):
        raise AbstractMethodError(self)


class PyArrowImpl(BaseImpl):
    def __init__(self):
        import_optional_dependency(
            "pyarrow", extra="pyarrow is required for parquet support."
        )
        import pyarrow.parquet

        # import utils to register the pyarrow extension types
        import pandas.core.arrays._arrow_utils  # noqa

        self.api = pyarrow

    def write(
        self,
        df: DataFrame,
        path: FilePathOrBuffer[AnyStr],
        compression: Optional[str] = "snappy",
        index: Optional[bool] = None,
        partition_cols: Optional[List[str]] = None,
        **kwargs,
    ):
        self.validate_dataframe(df)

        from_pandas_kwargs: Dict[str, Any] = {"schema": kwargs.pop("schema", None)}
        if index is not None:
            from_pandas_kwargs["preserve_index"] = index

        table = self.api.Table.from_pandas(df, **from_pandas_kwargs)

        if is_fsspec_url(path) and "filesystem" not in kwargs:
            # make fsspec instance, which pyarrow will use to open paths
            import_optional_dependency("fsspec")
            import fsspec.core

            fs, path = fsspec.core.url_to_fs(path)
            kwargs["filesystem"] = fs
        else:
            path = _expand_user(path)
        if partition_cols is not None:
            # writes to multiple files under the given path
            self.api.parquet.write_to_dataset(
                table,
                path,
                compression=compression,
                partition_cols=partition_cols,
                **kwargs,
            )
        else:
            # write to single output file
            self.api.parquet.write_table(table, path, compression=compression, **kwargs)

    def read(self, path, columns=None, **kwargs):
        if is_fsspec_url(path) and "filesystem" not in kwargs:
            import_optional_dependency("fsspec")
            import fsspec.core

            fs, path = fsspec.core.url_to_fs(path)
            should_close = False
        else:
            fs = kwargs.pop("filesystem", None)
            should_close = False
            path = _expand_user(path)

        if not fs:
            path, _, _, should_close = get_filepath_or_buffer(path)

        kwargs["use_pandas_metadata"] = True
        result = self.api.parquet.read_table(
            path, columns=columns, filesystem=fs, **kwargs
        ).to_pandas()
        if should_close:
            path.close()

        return result


class FastParquetImpl(BaseImpl):
    def __init__(self):
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
        compression="snappy",
        index=None,
        partition_cols=None,
        **kwargs,
    ):
        self.validate_dataframe(df)
        # thriftpy/protocol/compact.py:339:
        # DeprecationWarning: tostring() is deprecated.
        # Use tobytes() instead.

        if "partition_on" in kwargs and partition_cols is not None:
            raise ValueError(
                "Cannot use both partition_on and "
                "partition_cols. Use partition_cols for partitioning data"
            )
        elif "partition_on" in kwargs:
            partition_cols = kwargs.pop("partition_on")

        if partition_cols is not None:
            kwargs["file_scheme"] = "hive"

        if is_fsspec_url(path):
            fsspec = import_optional_dependency("fsspec")

            # if filesystem is provided by fsspec, file must be opened in 'wb' mode.
            kwargs["open_with"] = lambda path, _: fsspec.open(path, "wb").open()
        else:
            path, _, _, _ = get_filepath_or_buffer(path)

        with catch_warnings(record=True):
            self.api.write(
                path,
                df,
                compression=compression,
                write_index=index,
                partition_on=partition_cols,
                **kwargs,
            )

    def read(self, path, columns=None, **kwargs):
        if is_fsspec_url(path):
            fsspec = import_optional_dependency("fsspec")

            open_with = lambda path, _: fsspec.open(path, "rb").open()
            parquet_file = self.api.ParquetFile(path, open_with=open_with)
        else:
            path, _, _, _ = get_filepath_or_buffer(path)
            parquet_file = self.api.ParquetFile(path)

        return parquet_file.to_pandas(columns=columns, **kwargs)


def to_parquet(
    df: DataFrame,
    path: FilePathOrBuffer[AnyStr],
    engine: str = "auto",
    compression: Optional[str] = "snappy",
    index: Optional[bool] = None,
    partition_cols: Optional[List[str]] = None,
    **kwargs,
):
    """
    Write a DataFrame to the parquet format.

    Parameters
    ----------
    df : DataFrame
    path : str or file-like object
        If a string, it will be used as Root Directory path
        when writing a partitioned dataset. By file-like object,
        we refer to objects with a write() method, such as a file handler
        (e.g. via builtin open function) or io.BytesIO. The engine
        fastparquet does not accept file-like objects.

        .. versionchanged:: 0.24.0

    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet library to use. If 'auto', then the option
        ``io.parquet.engine`` is used. The default ``io.parquet.engine``
        behavior is to try 'pyarrow', falling back to 'fastparquet' if
        'pyarrow' is unavailable.
    compression : {'snappy', 'gzip', 'brotli', None}, default 'snappy'
        Name of the compression to use. Use ``None`` for no compression.
    index : bool, default None
        If ``True``, include the dataframe's index(es) in the file output. If
        ``False``, they will not be written to the file.
        If ``None``, similar to ``True`` the dataframe's index(es)
        will be saved. However, instead of being saved as values,
        the RangeIndex will be stored as a range in the metadata so it
        doesn't require much space and is faster. Other indexes will
        be included as columns in the file output.

        .. versionadded:: 0.24.0

    partition_cols : str or list, optional, default None
        Column names by which to partition the dataset.
        Columns are partitioned in the order they are given.
        Must be None if path is not a string.

        .. versionadded:: 0.24.0

    kwargs
        Additional keyword arguments passed to the engine
    """
    if isinstance(partition_cols, str):
        partition_cols = [partition_cols]
    impl = get_engine(engine)
    return impl.write(
        df,
        path,
        compression=compression,
        index=index,
        partition_cols=partition_cols,
        **kwargs,
    )


def read_parquet(path, engine: str = "auto", columns=None, **kwargs):
    """
    Load a parquet object from the file path, returning a DataFrame.

    Parameters
    ----------
    path : str, path object or file-like object
        Any valid string path is acceptable. The string could be a URL. Valid
        URL schemes include http, ftp, s3, and file. For file URLs, a host is
        expected. A local file could be:
        ``file://localhost/path/to/table.parquet``.
        A file URL can also be a path to a directory that contains multiple
        partitioned parquet files. Both pyarrow and fastparquet support
        paths to directories as well as file URLs. A directory path could be:
        ``file://localhost/path/to/tables`` or ``s3://bucket/partition_dir``

        If you want to pass in a path object, pandas accepts any
        ``os.PathLike``.

        By file-like object, we refer to objects with a ``read()`` method,
        such as a file handler (e.g. via builtin ``open`` function)
        or ``StringIO``.
    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet library to use. If 'auto', then the option
        ``io.parquet.engine`` is used. The default ``io.parquet.engine``
        behavior is to try 'pyarrow', falling back to 'fastparquet' if
        'pyarrow' is unavailable.
    columns : list, default=None
        If not None, only these columns will be read from the file.
    **kwargs
        Any additional kwargs are passed to the engine.

    Returns
    -------
    DataFrame
    """
    impl = get_engine(engine)
    return impl.read(path, columns=columns, **kwargs)
