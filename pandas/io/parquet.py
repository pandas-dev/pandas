""" parquet compat """

from warnings import catch_warnings
from distutils.version import LooseVersion
from pandas import DataFrame, RangeIndex, Int64Index, get_option
from pandas.compat import range
from pandas.io.common import get_filepath_or_buffer


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

    if engine not in ['pyarrow', 'fastparquet']:
        raise ValueError("engine must be one of 'pyarrow', 'fastparquet'")

    if engine == 'pyarrow':
        return PyArrowImpl()
    elif engine == 'fastparquet':
        return FastParquetImpl()


class PyArrowImpl(object):

    def __init__(self):
        # since pandas is a dependency of pyarrow
        # we need to import on first use

        try:
            import pyarrow
            import pyarrow.parquet
        except ImportError:
            raise ImportError("pyarrow is required for parquet support\n\n"
                              "you can install via conda\n"
                              "conda install pyarrow -c conda-forge\n"
                              "\nor via pip\n"
                              "pip install -U pyarrow\n")

        if LooseVersion(pyarrow.__version__) < '0.4.1':
            raise ImportError("pyarrow >= 0.4.1 is required for parquet"
                              "support\n\n"
                              "you can install via conda\n"
                              "conda install pyarrow -c conda-forge\n"
                              "\nor via pip\n"
                              "pip install -U pyarrow\n")

        self._pyarrow_lt_050 = LooseVersion(pyarrow.__version__) < '0.5.0'
        self._pyarrow_lt_060 = LooseVersion(pyarrow.__version__) < '0.6.0'
        self.api = pyarrow

    def write(self, df, path, compression='snappy',
              coerce_timestamps='ms', **kwargs):
        path, _, _ = get_filepath_or_buffer(path)
        if self._pyarrow_lt_060:
            table = self.api.Table.from_pandas(df, timestamps_to_ms=True)
            self.api.parquet.write_table(
                table, path, compression=compression, **kwargs)

        else:
            table = self.api.Table.from_pandas(df)
            self.api.parquet.write_table(
                table, path, compression=compression,
                coerce_timestamps=coerce_timestamps, **kwargs)

    def read(self, path):
        path, _, _ = get_filepath_or_buffer(path)
        return self.api.parquet.read_table(path).to_pandas()


class FastParquetImpl(object):

    def __init__(self):
        # since pandas is a dependency of fastparquet
        # we need to import on first use

        try:
            import fastparquet
        except ImportError:
            raise ImportError("fastparquet is required for parquet support\n\n"
                              "you can install via conda\n"
                              "conda install fastparquet -c conda-forge\n"
                              "\nor via pip\n"
                              "pip install -U fastparquet")

        if LooseVersion(fastparquet.__version__) < '0.1.0':
            raise ImportError("fastparquet >= 0.1.0 is required for parquet "
                              "support\n\n"
                              "you can install via conda\n"
                              "conda install fastparquet -c conda-forge\n"
                              "\nor via pip\n"
                              "pip install -U fastparquet")

        self.api = fastparquet

    def write(self, df, path, compression='snappy', **kwargs):
        # thriftpy/protocol/compact.py:339:
        # DeprecationWarning: tostring() is deprecated.
        # Use tobytes() instead.
        path, _, _ = get_filepath_or_buffer(path)
        with catch_warnings(record=True):
            self.api.write(path, df,
                           compression=compression, **kwargs)

    def read(self, path):
        path, _, _ = get_filepath_or_buffer(path)
        return self.api.ParquetFile(path).to_pandas()


def to_parquet(df, path, engine='auto', compression='snappy', **kwargs):
    """
    Write a DataFrame to the parquet format.

    Parameters
    ----------
    df : DataFrame
    path : string
        File path
    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet reader library to use. If 'auto', then the option
        'io.parquet.engine' is used. If 'auto', then the first
        library to be installed is used.
    compression : str, optional, default 'snappy'
        compression method, includes {'gzip', 'snappy', 'brotli'}
    kwargs
        Additional keyword arguments passed to the engine
    """

    impl = get_engine(engine)

    if not isinstance(df, DataFrame):
        raise ValueError("to_parquet only support IO with DataFrames")

    valid_types = {'string', 'unicode'}

    # validate index
    # --------------

    # validate that we have only a default index
    # raise on anything else as we don't serialize the index

    if not isinstance(df.index, Int64Index):
        raise ValueError("parquet does not support serializing {} "
                         "for the index; you can .reset_index()"
                         "to make the index into column(s)".format(
                             type(df.index)))

    if not df.index.equals(RangeIndex.from_range(range(len(df)))):
        raise ValueError("parquet does not support serializing a "
                         "non-default index for the index; you "
                         "can .reset_index() to make the index "
                         "into column(s)")

    if df.index.name is not None:
        raise ValueError("parquet does not serialize index meta-data on a "
                         "default index")

    # validate columns
    # ----------------

    # must have value column names (strings only)
    if df.columns.inferred_type not in valid_types:
        raise ValueError("parquet must have string column names")

    return impl.write(df, path, compression=compression)


def read_parquet(path, engine='auto', **kwargs):
    """
    Load a parquet object from the file path, returning a DataFrame.

    .. versionadded 0.21.0

    Parameters
    ----------
    path : string
        File path
    engine : {'auto', 'pyarrow', 'fastparquet'}, default 'auto'
        Parquet reader library to use. If 'auto', then the option
        'io.parquet.engine' is used. If 'auto', then the first
        library to be installed is used.
    kwargs are passed to the engine

    Returns
    -------
    DataFrame

    """

    impl = get_engine(engine)
    return impl.read(path)
