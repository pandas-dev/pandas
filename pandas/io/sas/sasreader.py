"""
Read SAS sas7bdat or xport files.
"""

from __future__ import annotations

from abc import (
    ABC,
    abstractmethod,
)
from collections.abc import Iterator
from typing import (
    TYPE_CHECKING,
    Self,
    overload,
)

from pandas.util._decorators import set_module

from pandas.io.common import stringify_path

if TYPE_CHECKING:
    from collections.abc import Hashable
    from types import TracebackType

    from pandas._typing import (
        CompressionOptions,
        FilePath,
        ReadBuffer,
    )

    from pandas import DataFrame


@set_module("pandas.api.typing")
class SASReader(Iterator["DataFrame"], ABC):
    """
    Abstract class for XportReader and SAS7BDATReader.
    """

    @abstractmethod
    def read(self, nrows: int | None = None) -> DataFrame: ...

    @abstractmethod
    def close(self) -> None: ...

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()


@overload
def read_sas(
    filepath_or_buffer: FilePath | ReadBuffer[bytes],
    *,
    format: str | None = ...,
    index: Hashable | None = ...,
    encoding: str | None = ...,
    chunksize: int = ...,
    iterator: bool = ...,
    compression: CompressionOptions = ...,
) -> SASReader: ...


@overload
def read_sas(
    filepath_or_buffer: FilePath | ReadBuffer[bytes],
    *,
    format: str | None = ...,
    index: Hashable | None = ...,
    encoding: str | None = ...,
    chunksize: None = ...,
    iterator: bool = ...,
    compression: CompressionOptions = ...,
) -> DataFrame | SASReader: ...


@set_module("pandas")
def read_sas(
    filepath_or_buffer: FilePath | ReadBuffer[bytes],
    *,
    format: str | None = None,
    index: Hashable | None = None,
    encoding: str | None = None,
    chunksize: int | None = None,
    iterator: bool = False,
    compression: CompressionOptions = "infer",
) -> DataFrame | SASReader:
    """
    Read SAS files stored as either XPORT or SAS7BDAT format files.

    Parameters
    ----------
    filepath_or_buffer : str, path object, or file-like object
        String, path object (implementing ``os.PathLike[str]``), or file-like
        object implementing a binary ``read()`` function. The string could be
        a URL. Valid URL schemes include http, ftp, s3, and file. For file
        URLs, a host is expected. A local file could be:
        ``file://localhost/path/to/table.sas7bdat``.
    format : str {{'xport', 'sas7bdat'}} or None
        If None, file format is inferred from file extension. If 'xport' or
        'sas7bdat', uses the corresponding format.
    index : identifier of index column, defaults to None
        Identifier of column that should be used as index of the DataFrame.
    encoding : str, default is None
        Encoding for text data.  If None, text data are stored as raw bytes.
    chunksize : int
        Read file `chunksize` lines at a time, returns iterator.
    iterator : bool, defaults to False
        If True, returns an iterator for reading the file incrementally.
    compression : str or dict, default 'infer'
        For on-the-fly decompression of on-disk data. If 'infer' and
        'filepath_or_buffer' is path-like, then detect compression from the
        following extensions: '.gz', '.bz2', '.zip', '.xz', '.zst', '.tar',
        '.tar.gz', '.tar.xz' or '.tar.bz2' (otherwise no compression).
        Set to ``None`` for no decompression.
        Can also be a dict with key ``'method'`` set to one of {``'zip'``,
        ``'gzip'``, ``'bz2'``, ``'zstd'``, ``'xz'``, ``'tar'``} and other
        key-value pairs are forwarded to ``zipfile.ZipFile``,
        ``gzip.GzipFile``, ``bz2.BZ2File``, ``zstandard.ZstdCompressor``,
        ``lzma.LZMAFile`` or ``tarfile.TarFile``, respectively.
        As an example, the following could be passed for faster compression
        and to create a reproducible gzip archive:
        ``compression={'method': 'gzip', 'compresslevel': 1, 'mtime': 1}``.

    Returns
    -------
    DataFrame, SAS7BDATReader, or XportReader
        DataFrame if iterator=False and chunksize=None, else SAS7BDATReader
        or XportReader, file format is inferred from file extension.

    See Also
    --------
    read_csv : Read a comma-separated values (csv) file into a DataFrame.
    read_excel : Read an Excel file into a pandas DataFrame.
    read_spss : Read an SPSS file into a pandas DataFrame.
    read_orc : Load an ORC object into a pandas DataFrame.
    read_feather : Load a feather-format object into a pandas DataFrame.

    Examples
    --------
    >>> df = pd.read_sas("sas_data.sas7bdat")  # doctest: +SKIP
    """
    if format is None:
        buffer_error_msg = (
            "If this is a buffer object rather "
            "than a string name, you must specify a format string"
        )
        filepath_or_buffer = stringify_path(filepath_or_buffer)
        if not isinstance(filepath_or_buffer, str):
            raise ValueError(buffer_error_msg)
        fname = filepath_or_buffer.lower()
        if ".xpt" in fname:
            format = "xport"
        elif ".sas7bdat" in fname:
            format = "sas7bdat"
        else:
            raise ValueError(
                f"unable to infer format of SAS file from filename: {fname!r}"
            )

    reader: SASReader
    if format.lower() == "xport":
        from pandas.io.sas.sas_xport import XportReader

        reader = XportReader(
            filepath_or_buffer,
            index=index,
            encoding=encoding,
            chunksize=chunksize,
            compression=compression,
        )
    elif format.lower() == "sas7bdat":
        from pandas.io.sas.sas7bdat import SAS7BDATReader

        reader = SAS7BDATReader(
            filepath_or_buffer,
            index=index,
            encoding=encoding,
            chunksize=chunksize,
            compression=compression,
        )
    else:
        raise ValueError("unknown SAS format")

    if iterator or chunksize:
        return reader

    with reader:
        return reader.read()
