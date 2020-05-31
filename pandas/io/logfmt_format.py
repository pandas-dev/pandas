"""logfmt format support"""

from collections import abc
from io import StringIO
from itertools import islice
from typing import Optional, Dict, Iterable, Generator, Union

from pandas._typing import FilePathOrBuffer, Dtype
from pandas.compat._optional import import_optional_dependency

from pandas import DataFrame
from pandas.core.indexes.api import RangeIndex
from pandas.core.reshape.concat import concat

from pandas.io.common import get_filepath_or_buffer, get_handle, infer_compression
from pandas.io.parsers import _validate_integer


def read_logfmt(
    filepath_or_buffer: FilePathOrBuffer,
    dtype: Optional[Dtype] = None,
    chunksize: Optional[int] = None,
    encoding: Optional[str] = None,
    compression: Optional[str] = "infer",
):
    """
    Load a logfmt_ from the file path.

    .. _logfmt: https://www.brandur.org/logfmt

    Parameters
    ----------
    filepath_or_buffer : path object or file-like object
        Any valid string path is acceptable. The string could be a URL.

        If you want to pass in a path object, pandas accepts any
        ``os.PathLike``.

        By file-like object, we refer to objects with a ``read()`` method,
        such as a file handler (e.g. via builtin ``open`` function)
        or ``StringIO``.

    dtype : bool or dict, default None
        If True, infer dtypes; if a dict of column to dtype, then use those;
        if False, then don't infer dtypes at all, applies only to the data.

    chunksize : int, optional
        Number of lines to be read per iteration.
        If None, read the whole file.

    encoding : str, optional
        Encoding to be used by the parser

    compression : str, optional
        Compression to be used by the parser

    Returns
    -------
    Data frame or `LogfmtReader`, if `chunksize` is specified.
    """

    compression = infer_compression(filepath_or_buffer, compression)
    filepath_or_buffer, _, _, should_close = get_filepath_or_buffer(
        filepath_or_buffer, encoding=encoding, compression=compression
    )

    logfmt_reader = LogfmtReader(
        filepath_or_buffer,
        dtype=dtype,
        chunksize=chunksize,
        encoding=encoding,
        compression=compression,
    )

    if chunksize:
        return logfmt_reader

    result = logfmt_reader.read()
    if should_close:
        result.close()

    return result


class LogfmtReader(abc.Iterator):
    """
    LogfmtReader provides an interface for reading in a logfmt file.
    """

    def __init__(
        self,
        filepath_or_buffer: FilePathOrBuffer,
        dtype: Optional[Dtype],
        chunksize: Optional[int],
        encoding: Optional[str],
        compression: Optional[str],
    ) -> None:
        self.chunksize = chunksize
        self.dtype = dtype
        self.encoding = encoding
        self.compression = compression
        self.nrows_seen = 0
        self.should_close = False

        if self.chunksize is not None:
            self.chunksize = _validate_integer("chunksize", self.chunksize, 1)

        if isinstance(filepath_or_buffer, str) or self.compression is not None:
            self.data, _ = get_handle(
                filepath_or_buffer,
                "r",
                encoding=self.encoding,
                compression=self.compression,
            )
            self.should_close = True
        else:
            self.data = filepath_or_buffer

    def _get_data_from_filepath(self, path: FilePathOrBuffer):
        return path

    def read(self) -> DataFrame:
        return concat(self)

    def close(self) -> None:
        """
        If we opened a stream earlier, we should close it.

        If an open stream or file was passed, we leave it open.
        """
        if self.should_close:
            try:
                self.data.close()
            except (IOError, AttributeError):
                pass

    @staticmethod
    def infer_types(lines: Iterable[Dict]) -> Generator[Dict, None, None]:
        """Infer types for parsed logfmt lines"""
        for line in lines:
            for key in line:
                try:
                    line[key] = int(line[key])
                except ValueError:
                    try:
                        line[key] = float(line[key])
                    except ValueError:
                        pass
            yield line

    def __next__(self):
        logfmt = import_optional_dependency("logfmt")

        lines = list(islice(self.data, self.chunksize))
        if lines:
            logfmt_lines = self.infer_types(logfmt.parse(StringIO("\n".join(lines))))
            obj = DataFrame(logfmt_lines, dtype=self.dtype)

            # Make sure that the returned objects have the right index.
            obj.index = RangeIndex(self.nrows_seen, self.nrows_seen + len(obj))
            self.nrows_seen += len(obj)

            return obj

        self.close()
        raise StopIteration
