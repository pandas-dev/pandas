from __future__ import annotations

from contextlib import contextmanager
import os
import sys
from typing import (
    IO,
    TYPE_CHECKING,
)

from pandas.compat import CHAINED_WARNING_DISABLED
from pandas.errors import ChainedAssignmentError

from pandas.io.common import get_handle

if TYPE_CHECKING:
    from collections.abc import Generator

    from pandas._typing import (
        BaseBuffer,
        CompressionOptions,
        FilePath,
    )


@contextmanager
def decompress_file(
    path: FilePath | BaseBuffer, compression: CompressionOptions
) -> Generator[IO[bytes]]:
    """
    Open a compressed file and return a file object.

    Parameters
    ----------
    path : str
        The path where the file is read from.

    compression : {'gzip', 'bz2', 'zip', 'xz', 'zstd', None}
        Name of the decompression to use

    Returns
    -------
    file object
    """
    with get_handle(path, "rb", compression=compression, is_text=False) as handle:
        yield handle.handle


@contextmanager
def set_timezone(tz: str) -> Generator[None]:
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
    >>> tzlocal().tzname(datetime(2021, 1, 1))  # doctest: +SKIP
    'IST'

    >>> with set_timezone("US/Eastern"):
    ...     tzlocal().tzname(datetime(2021, 1, 1))
    'EST'
    """
    import time

    def setTZ(tz) -> None:
        if hasattr(time, "tzset"):
            if tz is None:
                try:
                    del os.environ["TZ"]
                except KeyError:
                    pass
            else:
                os.environ["TZ"] = tz
                # Next line allows typing checks to pass on Windows
                if sys.platform != "win32":
                    time.tzset()

    orig_tz = os.environ.get("TZ")
    setTZ(tz)
    try:
        yield
    finally:
        setTZ(orig_tz)


@contextmanager
def with_csv_dialect(name: str, **kwargs) -> Generator[None]:
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
    try:
        yield
    finally:
        csv.unregister_dialect(name)


def raises_chained_assignment_error(extra_warnings=(), extra_match=()):
    from pandas._testing import assert_produces_warning

    if CHAINED_WARNING_DISABLED:
        if not extra_warnings:
            from contextlib import nullcontext

            return nullcontext()
        else:
            return assert_produces_warning(
                extra_warnings,
                match=extra_match,
            )
    else:
        warning = ChainedAssignmentError
        match = (
            "A value is being set on a copy of a DataFrame or Series "
            "through chained assignment"
        )
        if extra_warnings:
            warning = (warning, *extra_warnings)  # type: ignore[assignment]
        return assert_produces_warning(
            warning,
            match=(match, *extra_match),
        )
