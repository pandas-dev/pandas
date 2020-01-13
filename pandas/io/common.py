"""Common IO api utilities"""

import bz2
from collections import abc
import gzip
from io import BufferedIOBase, BytesIO
import mmap
import os
import pathlib
from typing import IO, Any, AnyStr, Dict, List, Mapping, Optional, Tuple, Union
from urllib.parse import (  # noqa
    urlencode,
    urljoin,
    urlparse as parse_url,
    uses_netloc,
    uses_params,
    uses_relative,
)
import zipfile

from pandas._typing import FilePathOrBuffer
from pandas.compat import _get_lzma_file, _import_lzma
from pandas.errors import (  # noqa
    AbstractMethodError,
    DtypeWarning,
    EmptyDataError,
    ParserError,
    ParserWarning,
)

from pandas.core.dtypes.common import is_file_like

lzma = _import_lzma()


_VALID_URLS = set(uses_relative + uses_netloc + uses_params)
_VALID_URLS.discard("")


def is_url(url) -> bool:
    """
    Check to see if a URL has a valid protocol.

    Parameters
    ----------
    url : str or unicode

    Returns
    -------
    isurl : bool
        If `url` has a valid protocol return True otherwise False.
    """
    if not isinstance(url, str):
        return False
    return parse_url(url).scheme in _VALID_URLS


def _expand_user(
    filepath_or_buffer: FilePathOrBuffer[AnyStr],
) -> FilePathOrBuffer[AnyStr]:
    """Return the argument with an initial component of ~ or ~user
       replaced by that user's home directory.

    Parameters
    ----------
    filepath_or_buffer : object to be converted if possible

    Returns
    -------
    expanded_filepath_or_buffer : an expanded filepath or the
                                  input if not expandable
    """
    if isinstance(filepath_or_buffer, str):
        return os.path.expanduser(filepath_or_buffer)
    return filepath_or_buffer


def validate_header_arg(header) -> None:
    if isinstance(header, bool):
        raise TypeError(
            "Passing a bool to header is invalid. "
            "Use header=None for no header or "
            "header=int or list-like of ints to specify "
            "the row(s) making up the column names"
        )


def stringify_path(
    filepath_or_buffer: FilePathOrBuffer[AnyStr],
) -> FilePathOrBuffer[AnyStr]:
    """Attempt to convert a path-like object to a string.

    Parameters
    ----------
    filepath_or_buffer : object to be converted

    Returns
    -------
    str_filepath_or_buffer : maybe a string version of the object

    Notes
    -----
    Objects supporting the fspath protocol (python 3.6+) are coerced
    according to its __fspath__ method.

    For backwards compatibility with older pythons, pathlib.Path and
    py.path objects are specially coerced.

    Any other object is passed through unchanged, which includes bytes,
    strings, buffers, or anything else that's not even path-like.
    """
    if hasattr(filepath_or_buffer, "__fspath__"):
        # https://github.com/python/mypy/issues/1424
        return filepath_or_buffer.__fspath__()  # type: ignore
    elif isinstance(filepath_or_buffer, pathlib.Path):
        return str(filepath_or_buffer)
    return _expand_user(filepath_or_buffer)


def is_s3_url(url) -> bool:
    """Check for an s3, s3n, or s3a url"""
    if not isinstance(url, str):
        return False
    return parse_url(url).scheme in ["s3", "s3n", "s3a"]


def is_gcs_url(url) -> bool:
    """Check for a gcs url"""
    if not isinstance(url, str):
        return False
    return parse_url(url).scheme in ["gcs", "gs"]


def urlopen(*args, **kwargs):
    """
    Lazy-import wrapper for stdlib urlopen, as that imports a big chunk of
    the stdlib.
    """
    import urllib.request

    return urllib.request.urlopen(*args, **kwargs)


def get_filepath_or_buffer(
    filepath_or_buffer: FilePathOrBuffer,
    encoding: Optional[str] = None,
    compression: Optional[str] = None,
    mode: Optional[str] = None,
):
    """
    If the filepath_or_buffer is a url, translate and return the buffer.
    Otherwise passthrough.

    Parameters
    ----------
    filepath_or_buffer : a url, filepath (str, py.path.local or pathlib.Path),
                         or buffer
    compression : {{'gzip', 'bz2', 'zip', 'xz', None}}, optional
    encoding : the encoding to use to decode bytes, default is 'utf-8'
    mode : str, optional

    Returns
    -------
    tuple of ({a filepath_ or buffer or S3File instance},
              encoding, str,
              compression, str,
              should_close, bool)
    """
    filepath_or_buffer = stringify_path(filepath_or_buffer)

    if isinstance(filepath_or_buffer, str) and is_url(filepath_or_buffer):
        req = urlopen(filepath_or_buffer)
        content_encoding = req.headers.get("Content-Encoding", None)
        if content_encoding == "gzip":
            # Override compression based on Content-Encoding header
            compression = "gzip"
        reader = BytesIO(req.read())
        req.close()
        return reader, encoding, compression, True

    if is_s3_url(filepath_or_buffer):
        from pandas.io import s3

        return s3.get_filepath_or_buffer(
            filepath_or_buffer, encoding=encoding, compression=compression, mode=mode
        )

    if is_gcs_url(filepath_or_buffer):
        from pandas.io import gcs

        return gcs.get_filepath_or_buffer(
            filepath_or_buffer, encoding=encoding, compression=compression, mode=mode
        )

    if isinstance(filepath_or_buffer, (str, bytes, mmap.mmap)):
        return _expand_user(filepath_or_buffer), None, compression, False

    if not is_file_like(filepath_or_buffer):
        msg = f"Invalid file path or buffer object type: {type(filepath_or_buffer)}"
        raise ValueError(msg)

    return filepath_or_buffer, None, compression, False


def file_path_to_url(path: str) -> str:
    """
    converts an absolute native path to a FILE URL.

    Parameters
    ----------
    path : a path in native format

    Returns
    -------
    a valid FILE URL
    """
    # lazify expensive import (~30ms)
    from urllib.request import pathname2url

    return urljoin("file:", pathname2url(path))


_compression_to_extension = {"gzip": ".gz", "bz2": ".bz2", "zip": ".zip", "xz": ".xz"}


def get_compression_method(
    compression: Optional[Union[str, Mapping[str, str]]]
) -> Tuple[Optional[str], Dict[str, str]]:
    """
    Simplifies a compression argument to a compression method string and
    a mapping containing additional arguments.

    Parameters
    ----------
    compression : str or mapping
        If string, specifies the compression method. If mapping, value at key
        'method' specifies compression method.

    Returns
    -------
    tuple of ({compression method}, Optional[str]
              {compression arguments}, Dict[str, str])

    Raises
    ------
    ValueError on mapping missing 'method' key
    """
    if isinstance(compression, Mapping):
        compression_args = dict(compression)
        try:
            compression = compression_args.pop("method")
        except KeyError:
            raise ValueError("If mapping, compression must have key 'method'")
    else:
        compression_args = {}
    return compression, compression_args


def infer_compression(
    filepath_or_buffer: FilePathOrBuffer, compression: Optional[str]
) -> Optional[str]:
    """
    Get the compression method for filepath_or_buffer. If compression='infer',
    the inferred compression method is returned. Otherwise, the input
    compression method is returned unchanged, unless it's invalid, in which
    case an error is raised.

    Parameters
    ----------
    filepath_or_buffer : str or file handle
        File path or object.
    compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}
        If 'infer' and `filepath_or_buffer` is path-like, then detect
        compression from the following extensions: '.gz', '.bz2', '.zip',
        or '.xz' (otherwise no compression).

    Returns
    -------
    string or None

    Raises
    ------
    ValueError on invalid compression specified.
    """

    # No compression has been explicitly specified
    if compression is None:
        return None

    # Infer compression
    if compression == "infer":
        # Convert all path types (e.g. pathlib.Path) to strings
        filepath_or_buffer = stringify_path(filepath_or_buffer)
        if not isinstance(filepath_or_buffer, str):
            # Cannot infer compression of a buffer, assume no compression
            return None

        # Infer compression from the filename/URL extension
        for compression, extension in _compression_to_extension.items():
            if filepath_or_buffer.endswith(extension):
                return compression
        return None

    # Compression has been specified. Check that it's valid
    if compression in _compression_to_extension:
        return compression

    msg = f"Unrecognized compression type: {compression}"
    valid = ["infer", None] + sorted(_compression_to_extension)
    msg += f"\nValid compression types are {valid}"
    raise ValueError(msg)


def get_handle(
    path_or_buf,
    mode: str,
    encoding=None,
    compression: Optional[Union[str, Mapping[str, Any]]] = None,
    memory_map: bool = False,
    is_text: bool = True,
):
    """
    Get file handle for given path/buffer and mode.

    Parameters
    ----------
    path_or_buf : str or file handle
        File path or object.
    mode : str
        Mode to open path_or_buf with.
    encoding : str or None
        Encoding to use.
    compression : str or dict, default None
        If string, specifies compression mode. If dict, value at key 'method'
        specifies compression mode. Compression mode must be one of {'infer',
        'gzip', 'bz2', 'zip', 'xz', None}. If compression mode is 'infer'
        and `filepath_or_buffer` is path-like, then detect compression from
        the following extensions: '.gz', '.bz2', '.zip', or '.xz' (otherwise
        no compression). If dict and compression mode is 'zip' or inferred as
        'zip', other entries passed as additional compression options.

        .. versionchanged:: 1.0.0

           May now be a dict with key 'method' as compression mode
           and other keys as compression options if compression
           mode is 'zip'.

    memory_map : boolean, default False
        See parsers._parser_params for more information.
    is_text : boolean, default True
        whether file/buffer is in text format (csv, json, etc.), or in binary
        mode (pickle, etc.).

    Returns
    -------
    f : file-like
        A file-like object.
    handles : list of file-like objects
        A list of file-like object that were opened in this function.
    """
    try:
        from s3fs import S3File

        need_text_wrapping = (BufferedIOBase, S3File)
    except ImportError:
        need_text_wrapping = BufferedIOBase  # type: ignore

    handles: List[IO] = list()
    f = path_or_buf

    # Convert pathlib.Path/py.path.local or string
    path_or_buf = stringify_path(path_or_buf)
    is_path = isinstance(path_or_buf, str)

    compression, compression_args = get_compression_method(compression)
    if is_path:
        compression = infer_compression(path_or_buf, compression)

    if compression:

        # GZ Compression
        if compression == "gzip":
            if is_path:
                f = gzip.open(path_or_buf, mode)
            else:
                f = gzip.GzipFile(fileobj=path_or_buf)

        # BZ Compression
        elif compression == "bz2":
            if is_path:
                f = bz2.BZ2File(path_or_buf, mode)
            else:
                f = bz2.BZ2File(path_or_buf)

        # ZIP Compression
        elif compression == "zip":
            zf = _BytesZipFile(path_or_buf, mode, **compression_args)
            # Ensure the container is closed as well.
            handles.append(zf)
            if zf.mode == "w":
                f = zf
            elif zf.mode == "r":
                zip_names = zf.namelist()
                if len(zip_names) == 1:
                    f = zf.open(zip_names.pop())
                elif len(zip_names) == 0:
                    raise ValueError(f"Zero files found in ZIP file {path_or_buf}")
                else:
                    raise ValueError(
                        "Multiple files found in ZIP file."
                        f" Only one file per ZIP: {zip_names}"
                    )

        # XZ Compression
        elif compression == "xz":
            f = _get_lzma_file(lzma)(path_or_buf, mode)

        # Unrecognized Compression
        else:
            msg = f"Unrecognized compression type: {compression}"
            raise ValueError(msg)

        handles.append(f)

    elif is_path:
        if encoding:
            # Encoding
            f = open(path_or_buf, mode, encoding=encoding, newline="")
        elif is_text:
            # No explicit encoding
            f = open(path_or_buf, mode, errors="replace", newline="")
        else:
            # Binary mode
            f = open(path_or_buf, mode)
        handles.append(f)

    # Convert BytesIO or file objects passed with an encoding
    if is_text and (compression or isinstance(f, need_text_wrapping)):
        from io import TextIOWrapper

        g = TextIOWrapper(f, encoding=encoding, newline="")
        if not isinstance(f, BufferedIOBase):
            handles.append(g)
        f = g

    if memory_map and hasattr(f, "fileno"):
        try:
            wrapped = _MMapWrapper(f)
            f.close()
            f = wrapped
        except Exception:
            # we catch any errors that may have occurred
            # because that is consistent with the lower-level
            # functionality of the C engine (pd.read_csv), so
            # leave the file handler as is then
            pass

    return f, handles


class _BytesZipFile(zipfile.ZipFile, BytesIO):  # type: ignore
    """
    Wrapper for standard library class ZipFile and allow the returned file-like
    handle to accept byte strings via `write` method.

    BytesIO provides attributes of file-like object and ZipFile.writestr writes
    bytes strings into a member of the archive.
    """

    # GH 17778
    def __init__(
        self,
        file: FilePathOrBuffer,
        mode: str,
        archive_name: Optional[str] = None,
        **kwargs,
    ):
        if mode in ["wb", "rb"]:
            mode = mode.replace("b", "")
        self.archive_name = archive_name
        super().__init__(file, mode, zipfile.ZIP_DEFLATED, **kwargs)

    def write(self, data):
        archive_name = self.filename
        if self.archive_name is not None:
            archive_name = self.archive_name
        super().writestr(archive_name, data)

    @property
    def closed(self):
        return self.fp is None


class _MMapWrapper(abc.Iterator):
    """
    Wrapper for the Python's mmap class so that it can be properly read in
    by Python's csv.reader class.

    Parameters
    ----------
    f : file object
        File object to be mapped onto memory. Must support the 'fileno'
        method or have an equivalent attribute

    """

    def __init__(self, f: IO):
        self.mmap = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)

    def __getattr__(self, name: str):
        return getattr(self.mmap, name)

    def __iter__(self) -> "_MMapWrapper":
        return self

    def __next__(self) -> str:
        newbytes = self.mmap.readline()

        # readline returns bytes, not str, but Python's CSV reader
        # expects str, so convert the output to str before continuing
        newline = newbytes.decode("utf-8")

        # mmap doesn't raise if reading past the allocated
        # data but instead returns an empty string, so raise
        # if that is returned
        if newline == "":
            raise StopIteration
        return newline
