""" parquet compat """

from distutils.version import LooseVersion
from warnings import catch_warnings

from pandas.errors import AbstractMethodError

from pandas import DataFrame, get_option

from pandas.io.common import get_filepath_or_buffer, is_s3_url


def get_engine(engine):
    """ return our implementation """

    if engine == 'auto':
        engine = get_option('io.parquet.engine')

    if engine == 'auto':
        # try engines in this order
        try:
            return PyArrowImpl()
        except ImportError:
            pass

        try:
            return FastParquetImpl()
        except ImportError:
            pass

        raise ImportError("Unable to find a usable engine; "
                          "tried using: 'pyarrow', 'fastparquet'.\n"
                          "pyarrow or fastparquet is required for parquet "
                          "support")

    if engine not in ['pyarrow', 'fastparquet']:
        raise ValueError("engine must be one of 'pyarrow', 'fastparquet'")

    if engine == 'pyarrow':
        return PyArrowImpl()
    elif engine == 'fastparquet':
        return FastParquetImpl()


class BaseImpl(object):

    api = None  # module

    @staticmethod
    def validate_dataframe(df):

        if not isinstance(df, DataFrame):
            raise ValueError("to_parquet only supports IO with DataFrames")

        # must have value column names (strings only)
        if df.columns.inferred_type not in {'string', 'unicode'}:
            raise ValueError("parquet must have string column names")

        # index level names must be strings
        valid_names = all(
            isinstance(name, str)
            for name in df.index.names
            if name is not None
        )
        if not valid_names:
            raise ValueError("Index level names must be strings")

    def write(self, df, path, compression, **kwargs):
        raise AbstractMethodError(self)

    def read(self, path, columns=None, **kwargs):
        raise AbstractMethodError(self)


class PyArrowImpl(BaseImpl):

    def __init__(self):
        # since pandas is a dependency of pyarrow
        # we need to import on first use
        try:
            import pyarrow
            import pyarrow.parquet
        except ImportError:
            raise ImportError(
                "pyarrow is required for parquet support\n\n"
                "you can install via conda\n"
                "conda install pyarrow -c conda-forge\n"
                "\nor via pip\n"
                "pip install -U pyarrow\n"
            )
        if LooseVersion(pyarrow.__version__) < '0.9.0':
            raise ImportError(
                "pyarrow >= 0.9.0 is required for parquet support\n\n"
                "you can install via conda\n"
                "conda install pyarrow -c conda-forge\n"
                "\nor via pip\n"
                "pip install -U pyarrow\n"
            )

        self.api = pyarrow

    def write(self, df, path, compression='snappy',
              coerce_timestamps='ms', index=None, partition_cols=None,
              **kwargs):
        self.validate_dataframe(df)
        path, _, _, _ = get_filepath_or_buffer(path, mode='wb')

        if index is None:
            from_pandas_kwargs = {}
        else:
            from_pandas_kwargs = {'preserve_index': index}
        table = self.api.Table.from_pandas(df, **from_pandas_kwargs)
        if partition_cols is not None:
            self.api.parquet.write_to_dataset(
                table, path, compression=compression,
                coerce_timestamps=coerce_timestamps,
                partition_cols=partition_cols, **kwargs)
        else:
            self.api.parquet.write_table(
                table, path, compression=compression,
                coerce_timestamps=coerce_timestamps, **kwargs)

    def read(self, path, columns=None, **kwargs):
        path, _, _, should_close = get_filepath_or_buffer(path)

        kwargs['use_pandas_metadata'] = True
        result = self.api.parquet.read_table(path, columns=columns,
                                             **kwargs).to_pandas()
        if should_close:
            try:
                path.close()
            except:  # noqa: flake8
                pass

        return result


class FastParquetImpl(BaseImpl):

    def __init__(self):
        # since pandas is a dependency of fastparquet
        # we need to import on first use
        try:
            import fastparquet
        except ImportError:
            raise ImportError(
                "fastparquet is required for parquet support\n\n"
                "you can install via conda\n"
                "conda install fastparquet -c conda-forge\n"
                "\nor via pip\n"
                "pip install -U fastparquet"
            )
        if LooseVersion(fastparquet.__version__) < '0.2.1':
            raise ImportError(
                "fastparquet >= 0.2.1 is required for parquet "
                "support\n\n"
                "you can install via conda\n"
                "conda install fastparquet -c conda-forge\n"
                "\nor via pip\n"
                "pip install -U fastparquet"
            )
        self.api = fastparquet

    def write(self, df, path, compression='snappy', index=None,
              partition_cols=None, **kwargs):
        self.validate_dataframe(df)
        # thriftpy/protocol/compact.py:339:
        # DeprecationWarning: tostring() is deprecated.
        # Use tobytes() instead.

        if 'partition_on' in kwargs and partition_cols is not None:
            raise ValueError("Cannot use both partition_on and "
                             "partition_cols. Use partition_cols for "
                             "partitioning data")
        elif 'partition_on' in kwargs:
            partition_cols = kwargs.pop('partition_on')

        if partition_cols is not None:
            kwargs['file_scheme'] = 'hive'

        if is_s3_url(path):
            # path is s3:// so we need to open the s3file in 'wb' mode.
            # TODO: Support 'ab'

            path, _, _, _ = get_filepath_or_buffer(path, mode='wb')
            # And pass the opened s3file to the fastparquet internal impl.
            kwargs['open_with'] = lambda path, _: path
        else:
            path, _, _, _ = get_filepath_or_buffer(path)

        with catch_warnings(record=True):
            self.api.write(path, df, compression=compression,
                           write_index=index, partition_on=partition_cols,
                           **kwargs)

    def read(self, path, columns=None, **kwargs):
        if is_s3_url(path):
            # When path is s3:// an S3File is returned.
            # We need to retain the original path(str) while also
            # pass the S3File().open function to fsatparquet impl.
            s3, _, _, should_close = get_filepath_or_buffer(path)
            try:
                parquet_file = self.api.ParquetFile(path, open_with=s3.s3.open)
            finally:
                s3.close()
        else:
            path, _, _, _ = get_filepath_or_buffer(path)
            parquet_file = self.api.ParquetFile(path)

        return parquet_file.to_pandas(columns=columns, **kwargs)


def to_parquet(df, path, engine='auto', compression='snappy', index=None,
               partition_cols=None, **kwargs):
    """
    Write a DataFrame to the parquet format.

    Parameters
    ----------
    path : str
        File path or Root Directory path. Will be used as Root Directory path
        while writing a partitioned dataset.

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
        ``False``, they will not be written to the file. If ``None``, the
        engine's default behavior will be used.

        .. versionadded 0.24.0

    partition_cols : list, optional, default None
        Column names by which to partition the dataset
        Columns are partitioned in the order they are given

        .. versionadded:: 0.24.0

    kwargs
        Additional keyword arguments passed to the engine
    """
    impl = get_engine(engine)
    return impl.write(df, path, compression=compression, index=index,
                      partition_cols=partition_cols, **kwargs)


def read_parquet(path, engine='auto', columns=None, **kwargs):
    """
    Load a parquet object from the file path, returning a DataFrame.

    .. versionadded 0.21.0

    Parameters
    ----------
    path : string
        File path
    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet library to use. If 'auto', then the option
        ``io.parquet.engine`` is used. The default ``io.parquet.engine``
        behavior is to try 'pyarrow', falling back to 'fastparquet' if
        'pyarrow' is unavailable.
    columns : list, default=None
        If not None, only these columns will be read from the file.

        .. versionadded 0.21.1
    **kwargs
        Any additional kwargs are passed to the engine.

    Returns
    -------
    DataFrame
    """

    impl = get_engine(engine)
    return impl.read(path, columns=columns, **kwargs)
