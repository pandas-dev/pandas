import bz2
from contextlib import contextmanager
import gzip
import os
from shutil import rmtree
import tempfile
import zipfile

from pandas.compat import get_lzma_file, import_lzma

lzma = import_lzma()


@contextmanager
def decompress_file(path, compression):
    """
    Open a compressed file and return a file object.

    Parameters
    ----------
    path : str
        The path where the file is read from.

    compression : {'gzip', 'bz2', 'zip', 'xz', None}
        Name of the decompression to use

    Returns
    -------
    file object
    """
    if compression is None:
        f = open(path, "rb")
    elif compression == "gzip":
        # pandas\_testing.py:243: error: Incompatible types in assignment
        # (expression has type "IO[Any]", variable has type "BinaryIO")
        f = gzip.open(path, "rb")  # type: ignore[assignment]
    elif compression == "bz2":
        # pandas\_testing.py:245: error: Incompatible types in assignment
        # (expression has type "BZ2File", variable has type "BinaryIO")
        f = bz2.BZ2File(path, "rb")  # type: ignore[assignment]
    elif compression == "xz":
        f = get_lzma_file(lzma)(path, "rb")
    elif compression == "zip":
        zip_file = zipfile.ZipFile(path)
        zip_names = zip_file.namelist()
        if len(zip_names) == 1:
            # pandas\_testing.py:252: error: Incompatible types in assignment
            # (expression has type "IO[bytes]", variable has type "BinaryIO")
            f = zip_file.open(zip_names.pop())  # type: ignore[assignment]
        else:
            raise ValueError(f"ZIP file {path} error. Only one file per ZIP.")
    else:
        raise ValueError(f"Unrecognized compression type: {compression}")

    try:
        yield f
    finally:
        f.close()
        if compression == "zip":
            zip_file.close()


@contextmanager
def set_timezone(tz: str):
    """
    Context manager for temporarily setting a timezone.

    Parameters
    ----------
    tz : str
        A string representing a valid timezone.

    Examples
    --------
    >>> from datetime import datetime
    >>> from dateutil.tz import tzlocal
    >>> tzlocal().tzname(datetime.now())
    'IST'

    >>> with set_timezone('US/Eastern'):
    ...     tzlocal().tzname(datetime.now())
    ...
    'EDT'
    """
    import os
    import time

    def setTZ(tz):
        if tz is None:
            try:
                del os.environ["TZ"]
            except KeyError:
                pass
        else:
            os.environ["TZ"] = tz
            time.tzset()

    orig_tz = os.environ.get("TZ")
    setTZ(tz)
    try:
        yield
    finally:
        setTZ(orig_tz)


@contextmanager
def ensure_clean(filename=None, return_filelike=False, **kwargs):
    """
    Gets a temporary path and agrees to remove on close.

    Parameters
    ----------
    filename : str (optional)
        if None, creates a temporary file which is then removed when out of
        scope. if passed, creates temporary file with filename as ending.
    return_filelike : bool (default False)
        if True, returns a file-like which is *always* cleaned. Necessary for
        savefig and other functions which want to append extensions.
    **kwargs
        Additional keywords passed in for creating a temporary file.
        :meth:`tempFile.TemporaryFile` is used when `return_filelike` is ``True``.
        :meth:`tempfile.mkstemp` is used when `return_filelike` is ``False``.
        Note that the `filename` parameter will be passed in as the `suffix`
        argument to either function.

    See Also
    --------
    tempfile.TemporaryFile
    tempfile.mkstemp
    """
    filename = filename or ""
    fd = None

    kwargs["suffix"] = filename

    if return_filelike:
        f = tempfile.TemporaryFile(**kwargs)

        try:
            yield f
        finally:
            f.close()
    else:
        # Don't generate tempfile if using a path with directory specified.
        if len(os.path.dirname(filename)):
            raise ValueError("Can't pass a qualified name to ensure_clean()")

        try:
            fd, filename = tempfile.mkstemp(**kwargs)
        except UnicodeEncodeError:
            import pytest

            pytest.skip("no unicode file names on this system")

        try:
            yield filename
        finally:
            try:
                os.close(fd)
            except OSError:
                print(f"Couldn't close file descriptor: {fd} (file: {filename})")
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except OSError as e:
                print(f"Exception on removing file: {e}")


@contextmanager
def ensure_clean_dir():
    """
    Get a temporary directory path and agrees to remove on close.

    Yields
    ------
    Temporary directory path
    """
    directory_name = tempfile.mkdtemp(suffix="")
    try:
        yield directory_name
    finally:
        try:
            rmtree(directory_name)
        except OSError:
            pass


@contextmanager
def ensure_safe_environment_variables():
    """
    Get a context manager to safely set environment variables

    All changes will be undone on close, hence environment variables set
    within this contextmanager will neither persist nor change global state.
    """
    saved_environ = dict(os.environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(saved_environ)


@contextmanager
def with_csv_dialect(name, **kwargs):
    """
    Context manager to temporarily register a CSV dialect for parsing CSV.

    Parameters
    ----------
    name : str
        The name of the dialect.
    kwargs : mapping
        The parameters for the dialect.

    Raises
    ------
    ValueError : the name of the dialect conflicts with a builtin one.

    See Also
    --------
    csv : Python's CSV library.
    """
    import csv

    _BUILTIN_DIALECTS = {"excel", "excel-tab", "unix"}

    if name in _BUILTIN_DIALECTS:
        raise ValueError("Cannot override builtin dialect.")

    csv.register_dialect(name, **kwargs)
    yield
    csv.unregister_dialect(name)


@contextmanager
def use_numexpr(use, min_elements=None):
    from pandas.core.computation import expressions as expr

    if min_elements is None:
        min_elements = expr._MIN_ELEMENTS

    olduse = expr.USE_NUMEXPR
    oldmin = expr._MIN_ELEMENTS
    expr.set_use_numexpr(use)
    expr._MIN_ELEMENTS = min_elements
    yield
    expr._MIN_ELEMENTS = oldmin
    expr.set_use_numexpr(olduse)
