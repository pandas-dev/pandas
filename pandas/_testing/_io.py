from __future__ import annotations

import gzip
import io
import tarfile
from typing import (
    TYPE_CHECKING,
    Any,
)
import zipfile

from pandas.compat._optional import import_optional_dependency

import pandas as pd

if TYPE_CHECKING:
    from collections.abc import Callable
    from pathlib import Path

    from pandas import (
        DataFrame,
        Series,
    )

# ------------------------------------------------------------------
# File-IO


def round_trip_pickle(obj: Any, tmp_path: Path) -> DataFrame | Series:
    """
    Pickle an object and then read it again.

    Parameters
    ----------
    obj : any object
        The object to pickle and then re-read.
    path : str, path object or file-like object, default None
        The path where the pickled object is written and then read.

    Returns
    -------
    pandas object
        The original object that was pickled and then re-read.
    """
    pd.to_pickle(obj, tmp_path)
    return pd.read_pickle(tmp_path)


def round_trip_pathlib(writer, reader, tmp_path: Path):
    """
    Write an object to file specified by a pathlib.Path and read it back

    Parameters
    ----------
    writer : callable bound to pandas object
        IO writing function (e.g. DataFrame.to_csv )
    reader : callable
        IO reading function (e.g. pd.read_csv )
    path : str, default None
        The path where the object is written and then read.

    Returns
    -------
    pandas object
        The original object that was serialized and then re-read.
    """
    writer(tmp_path)
    obj = reader(tmp_path)
    return obj


def write_to_compressed(compression, path: str, data, dest: str = "test") -> None:
    """
    Write data to a compressed file.

    Parameters
    ----------
    compression : {'gzip', 'bz2', 'zip', 'xz', 'zstd'}
        The compression type to use.
    path : str
        The file path to write the data.
    data : str
        The data to write.
    dest : str, default "test"
        The destination file (for ZIP only)

    Raises
    ------
    ValueError : An invalid compression value was passed in.
    """
    args: tuple[Any, ...] = (data,)
    mode = "wb"
    method = "write"
    compress_method: Callable

    if compression == "zip":
        compress_method = zipfile.ZipFile
        mode = "w"
        args = (dest, data)
        method = "writestr"
    elif compression == "tar":
        compress_method = tarfile.TarFile
        mode = "w"
        file = tarfile.TarInfo(name=dest)
        bytes = io.BytesIO(data)
        file.size = len(data)
        args = (file, bytes)
        method = "addfile"
    elif compression == "gzip":
        compress_method = gzip.GzipFile
    elif compression == "bz2":
        import bz2

        compress_method = bz2.BZ2File
    elif compression == "zstd":
        compress_method = import_optional_dependency("zstandard").open
    elif compression == "xz":
        import lzma

        compress_method = lzma.LZMAFile
    else:
        raise ValueError(f"Unrecognized compression type: {compression}")

    # error: No overload variant of "ZipFile" matches argument types "str", "str"
    # error: No overload variant of "BZ2File" matches argument types "str", "str"
    # error: Argument "mode" to "TarFile" has incompatible type "str";
    #  expected "Literal['r', 'a', 'w', 'x']
    with compress_method(path, mode=mode) as f:  # type: ignore[call-overload, arg-type]
        getattr(f, method)(*args)
