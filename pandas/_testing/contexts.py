import bz2
from contextlib import contextmanager
import gzip
import os
from shutil import rmtree
import tempfile
import warnings
import zipfile

import numpy as np

from pandas.compat import _get_lzma_file, _import_lzma

from pandas.core.dtypes.common import is_list_like

lzma = _import_lzma()


@contextmanager
def assert_produces_warning(
    expected_warning=Warning,
    filter_level="always",
    clear=None,
    check_stacklevel=True,
    raise_on_extra_warnings=True,
):
    """
    Context manager for running code expected to either raise a specific
    warning, or not raise any warnings. Verifies that the code raises the
    expected warning, and that it does not raise any other unexpected
    warnings. It is basically a wrapper around ``warnings.catch_warnings``.

    Parameters
    ----------
    expected_warning : {Warning, False, None}, default Warning
        The type of Exception raised. ``exception.Warning`` is the base
        class for all warnings. To check that no warning is returned,
        specify ``False`` or ``None``.
    filter_level : str or None, default "always"
        Specifies whether warnings are ignored, displayed, or turned
        into errors.
        Valid values are:

        * "error" - turns matching warnings into exceptions
        * "ignore" - discard the warning
        * "always" - always emit a warning
        * "default" - print the warning the first time it is generated
          from each location
        * "module" - print the warning the first time it is generated
          from each module
        * "once" - print the warning the first time it is generated

    clear : str, default None
        If not ``None`` then remove any previously raised warnings from
        the ``__warningsregistry__`` to ensure that no warning messages are
        suppressed by this context manager. If ``None`` is specified,
        the ``__warningsregistry__`` keeps track of which warnings have been
        shown, and does not show them again.
    check_stacklevel : bool, default True
        If True, displays the line that called the function containing
        the warning to show were the function is called. Otherwise, the
        line that implements the function is displayed.
    raise_on_extra_warnings : bool, default True
        Whether extra warnings not of the type `expected_warning` should
        cause the test to fail.

    Examples
    --------
    >>> import warnings
    >>> with assert_produces_warning():
    ...     warnings.warn(UserWarning())
    ...
    >>> with assert_produces_warning(False):
    ...     warnings.warn(RuntimeWarning())
    ...
    Traceback (most recent call last):
        ...
    AssertionError: Caused unexpected warning(s): ['RuntimeWarning'].
    >>> with assert_produces_warning(UserWarning):
    ...     warnings.warn(RuntimeWarning())
    Traceback (most recent call last):
        ...
    AssertionError: Did not see expected warning of class 'UserWarning'.

    ..warn:: This is *not* thread-safe.
    """
    __tracebackhide__ = True

    with warnings.catch_warnings(record=True) as w:

        if clear is not None:
            # make sure that we are clearing these warnings
            # if they have happened before
            # to guarantee that we will catch them
            if not is_list_like(clear):
                clear = [clear]
            for m in clear:
                try:
                    m.__warningregistry__.clear()
                except AttributeError:
                    # module may not have __warningregistry__
                    pass

        saw_warning = False
        warnings.simplefilter(filter_level)
        yield w
        extra_warnings = []

        for actual_warning in w:
            if expected_warning and issubclass(
                actual_warning.category, expected_warning
            ):
                saw_warning = True

                if check_stacklevel and issubclass(
                    actual_warning.category, (FutureWarning, DeprecationWarning)
                ):
                    from inspect import getframeinfo, stack

                    caller = getframeinfo(stack()[2][0])
                    msg = (
                        "Warning not set with correct stacklevel. "
                        f"File where warning is raised: {actual_warning.filename} != "
                        f"{caller.filename}. Warning message: {actual_warning.message}"
                    )
                    assert actual_warning.filename == caller.filename, msg
            else:
                extra_warnings.append(
                    (
                        actual_warning.category.__name__,
                        actual_warning.message,
                        actual_warning.filename,
                        actual_warning.lineno,
                    )
                )
        if expected_warning:
            msg = (
                f"Did not see expected warning of class "
                f"{repr(expected_warning.__name__)}"
            )
            assert saw_warning, msg
        if raise_on_extra_warnings and extra_warnings:
            raise AssertionError(
                f"Caused unexpected warning(s): {repr(extra_warnings)}"
            )


class RNGContext:
    """
    Context manager to set the numpy random number generator speed. Returns
    to the original value upon exiting the context manager.

    Parameters
    ----------
    seed : int
        Seed for numpy.random.seed

    Examples
    --------
    with RNGContext(42):
        np.random.randn()
    """

    def __init__(self, seed):
        self.seed = seed

    def __enter__(self):

        self.start_state = np.random.get_state()
        np.random.seed(self.seed)

    def __exit__(self, exc_type, exc_value, traceback):

        np.random.set_state(self.start_state)


@contextmanager
def use_numexpr(use, min_elements=None):
    from pandas.core.computation import expressions as expr

    if min_elements is None:
        min_elements = expr._MIN_ELEMENTS

    olduse = expr._USE_NUMEXPR
    oldmin = expr._MIN_ELEMENTS
    expr.set_use_numexpr(use)
    expr._MIN_ELEMENTS = min_elements
    yield
    expr._MIN_ELEMENTS = oldmin
    expr.set_use_numexpr(olduse)


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
        f = gzip.open(path, "rb")
    elif compression == "bz2":
        f = bz2.BZ2File(path, "rb")
    elif compression == "xz":
        f = _get_lzma_file(lzma)(path, "rb")
    elif compression == "zip":
        zip_file = zipfile.ZipFile(path)
        zip_names = zip_file.namelist()
        if len(zip_names) == 1:
            f = zip_file.open(zip_names.pop())
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


# -----------------------------------------------------------------------------
# contextmanager to ensure the file cleanup


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
