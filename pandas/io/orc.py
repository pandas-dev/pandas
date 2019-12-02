""" orc compat """

import distutils
from typing import TYPE_CHECKING, List, Optional

from pandas.compat._optional import import_optional_dependency

from pandas._typing import FilePathOrBuffer

from pandas.io.common import get_filepath_or_buffer

if TYPE_CHECKING:
    from pandas import DataFrame


class PyArrowImpl:
    def __init__(self):
        pyarrow = import_optional_dependency(
            "pyarrow", extra="pyarrow is required for orc support."
        )

        # we require a newer version of pyarrow thaN we support for parquet
        import pyarrow

        if distutils.version.LooseVersion(pyarrow.__version__) < "0.13.0":
            raise ImportError("pyarrow must be >= 0.13.0 for read_orc")

        import pyarrow.orc

        self.api = pyarrow

    def read(
        self, path: FilePathOrBuffer, columns: Optional[List[str]] = None, **kwargs
    ) -> "DataFrame":
        path, _, _, _ = get_filepath_or_buffer(path)
        orc_file = self.api.orc.ORCFile(path)
        result = orc_file.read(columns=columns, **kwargs).to_pandas()
        return result


def read_orc(
    path: FilePathOrBuffer, columns: Optional[List[str]] = None, **kwargs,
):
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
    columns : list, default=None
        If not None, only these columns will be read from the file.
    **kwargs
        Any additional kwargs are passed to pyarrow.

    Returns
    -------
    DataFrame
    """

    impl = PyArrowImpl()
    return impl.read(path, columns=columns, **kwargs)
