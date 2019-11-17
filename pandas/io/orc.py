""" orc compat """

from pandas.compat._optional import import_optional_dependency
from pandas.errors import AbstractMethodError

from pandas import DataFrame, get_option

from pandas.io.common import get_filepath_or_buffer


def get_engine(engine):
    """ return our implementation """

    if engine == "auto":
        engine = get_option("io.orc.engine")

    if engine == "auto":
        # try engines in this order
        try:
            return PyArrowImpl()
        except ImportError:
            pass

        raise ImportError(
            "Unable to find a usable engine; "
            "tried using: 'pyarrow'.\n"
            "pyarrow is required for orc "
            "support"
        )

    if engine not in ["pyarrow"]:
        raise ValueError("engine must be 'pyarrow'")

    if engine == "pyarrow":
        return PyArrowImpl()


class BaseImpl:

    api = None  # module

    @staticmethod
    def validate_dataframe(df):

        if not isinstance(df, DataFrame):
            raise ValueError("to_orc only supports IO with DataFrames")

        # must have value column names (strings only)
        if df.columns.inferred_type not in {"string", "unicode", "empty"}:
            raise ValueError("ORC must have string column names")

        # index level names must be strings
        valid_names = all(
            isinstance(name, str) for name in df.index.names if name is not None
        )
        if not valid_names:
            raise ValueError("Index level names must be strings")

    def write(self, df, path, compression, **kwargs):
        raise AbstractMethodError(self)

    def read(self, path, columns=None, **kwargs):
        raise AbstractMethodError(self)


class PyArrowImpl(BaseImpl):
    def __init__(self):
        pyarrow = import_optional_dependency(
            "pyarrow", extra="pyarrow is required for orc support."
        )
        import pyarrow.orc

        self.api = pyarrow

    def read(self, path, columns=None, **kwargs):
        path, _, _, _ = get_filepath_or_buffer(path)

        py_file = self.api.input_stream(path)
        orc_file = self.api.orc.ORCFile(py_file)

        result = orc_file.read(columns=columns, **kwargs).to_pandas()

        return result


def read_orc(path, engine="auto", columns=None, **kwargs):
    """
    Load an ORC object from the file path, returning a DataFrame.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    path : str, path object or file-like object
        Any valid string path is acceptable. The string could be a URL. Valid
        URL schemes include http, ftp, s3, and file. For file URLs, a host is
        expected. A local file could be:
        ``file://localhost/path/to/table.orc``.

        If you want to pass in a path object, pandas accepts any
        ``os.PathLike``.

        By file-like object, we refer to objects with a ``read()`` method,
        such as a file handler (e.g. via builtin ``open`` function)
        or ``StringIO``.
    engine : {'auto', 'pyarrow'}, default 'auto'
        ORC library to use. If 'auto', then the option ``io.orc.engine`` is
        used. The default ``io.orc.engine`` behavior is to try 'pyarrow'.
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
