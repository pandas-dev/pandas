"""
This module provides decorator functions which can be applied to test objects
in order to skip those objects when certain conditions occur. A sample use case
is to detect if the platform is missing ``matplotlib``. If so, any test objects
which require ``matplotlib`` and decorated with ``@td.skip_if_no_mpl`` will be
skipped by ``pytest`` during the execution of the test suite.

To illustrate, after importing this module:

import pandas.util._test_decorators as td

The decorators can be applied to classes:

@td.skip_if_some_reason
class Foo:
    ...

Or individual functions:

@td.skip_if_some_reason
def test_foo():
    ...

For more information, refer to the ``pytest`` documentation on ``skipif``.
"""
from distutils.version import LooseVersion
from functools import wraps
import locale
from typing import Callable, Optional

import numpy as np
import pytest

from pandas.compat import is_platform_32bit, is_platform_windows
from pandas.compat._optional import import_optional_dependency
from pandas.compat.numpy import _np_version

from pandas.core.computation.expressions import _NUMEXPR_INSTALLED, _USE_NUMEXPR


def safe_import(mod_name: str, min_version: Optional[str] = None):
    """
    Parameters:
    -----------
    mod_name : str
        Name of the module to be imported
    min_version : str, default None
        Minimum required version of the specified mod_name

    Returns:
    --------
    object
        The imported module if successful, or False
    """
    try:
        mod = __import__(mod_name)
    except ImportError:
        return False

    if not min_version:
        return mod
    else:
        import sys

        try:
            version = getattr(sys.modules[mod_name], "__version__")
        except AttributeError:
            # xlrd uses a capitalized attribute name
            version = getattr(sys.modules[mod_name], "__VERSION__")
        if version:
            from distutils.version import LooseVersion

            if LooseVersion(version) >= LooseVersion(min_version):
                return mod

    return False


# TODO:
# remove when gh-24839 is fixed.
# this affects numpy 1.16 and pytables 3.4.4
tables = safe_import("tables")
xfail_non_writeable = pytest.mark.xfail(
    tables
    and LooseVersion(np.__version__) >= LooseVersion("1.16")
    and LooseVersion(tables.__version__) < LooseVersion("3.5.1"),
    reason=(
        "gh-25511, gh-24839. pytables needs a "
        "release beyond 3.4.4 to support numpy 1.16.x"
    ),
)


def _skip_if_no_mpl():
    mod = safe_import("matplotlib")
    if mod:
        mod.use("Agg", warn=True)
    else:
        return True


def _skip_if_has_locale():
    lang, _ = locale.getlocale()
    if lang is not None:
        return True


def _skip_if_not_us_locale():
    lang, _ = locale.getlocale()
    if lang != "en_US":
        return True


def _skip_if_no_scipy() -> bool:
    return not (
        safe_import("scipy.stats")
        and safe_import("scipy.sparse")
        and safe_import("scipy.interpolate")
        and safe_import("scipy.signal")
    )


def skip_if_installed(package: str) -> Callable:
    """
    Skip a test if a package is installed.

    Parameters
    ----------
    package : str
        The name of the package.
    """
    return pytest.mark.skipif(
        safe_import(package), reason=f"Skipping because {package} is installed."
    )


def skip_if_no(package: str, min_version: Optional[str] = None) -> Callable:
    """
    Generic function to help skip tests when required packages are not
    present on the testing system.

    This function returns a pytest mark with a skip condition that will be
    evaluated during test collection. An attempt will be made to import the
    specified ``package`` and optionally ensure it meets the ``min_version``

    The mark can be used as either a decorator for a test function or to be
    applied to parameters in pytest.mark.parametrize calls or parametrized
    fixtures.

    If the import and version check are unsuccessful, then the test function
    (or test case when used in conjunction with parametrization) will be
    skipped.

    Parameters
    ----------
    package: str
        The name of the required package.
    min_version: str or None, default None
        Optional minimum version of the package.

    Returns
    -------
    _pytest.mark.structures.MarkDecorator
        a pytest.mark.skipif to use as either a test decorator or a
        parametrization mark.
    """
    msg = f"Could not import '{package}'"
    if min_version:
        msg += f" satisfying a min_version of {min_version}"
    return pytest.mark.skipif(
        not safe_import(package, min_version=min_version), reason=msg
    )


skip_if_no_mpl = pytest.mark.skipif(
    _skip_if_no_mpl(), reason="Missing matplotlib dependency"
)
skip_if_mpl = pytest.mark.skipif(not _skip_if_no_mpl(), reason="matplotlib is present")
skip_if_32bit = pytest.mark.skipif(is_platform_32bit(), reason="skipping for 32 bit")
skip_if_windows = pytest.mark.skipif(is_platform_windows(), reason="Running on Windows")
skip_if_windows_python_3 = pytest.mark.skipif(
    is_platform_windows(), reason="not used on win32"
)
skip_if_has_locale = pytest.mark.skipif(
    _skip_if_has_locale(), reason=f"Specific locale is set {locale.getlocale()[0]}",
)
skip_if_not_us_locale = pytest.mark.skipif(
    _skip_if_not_us_locale(), reason=f"Specific locale is set {locale.getlocale()[0]}",
)
skip_if_no_scipy = pytest.mark.skipif(
    _skip_if_no_scipy(), reason="Missing SciPy requirement"
)
skip_if_no_ne = pytest.mark.skipif(
    not _USE_NUMEXPR,
    reason=f"numexpr enabled->{_USE_NUMEXPR}, installed->{_NUMEXPR_INSTALLED}",
)


def skip_if_np_lt(
    ver_str: str, reason: Optional[str] = None, *args, **kwds
) -> Callable:
    if reason is None:
        reason = f"NumPy {ver_str} or greater required"
    return pytest.mark.skipif(
        _np_version < LooseVersion(ver_str), reason=reason, *args, **kwds
    )


def parametrize_fixture_doc(*args):
    """
    Intended for use as a decorator for parametrized fixture,
    this function will wrap the decorated function with a pytest
    ``parametrize_fixture_doc`` mark. That mark will format
    initial fixture docstring by replacing placeholders {0}, {1} etc
    with parameters passed as arguments.

    Parameters
    ----------
    args: iterable
        Positional arguments for docstring.

    Returns
    -------
    function
        The decorated function wrapped within a pytest
        ``parametrize_fixture_doc`` mark
    """

    def documented_fixture(fixture):
        fixture.__doc__ = fixture.__doc__.format(*args)
        return fixture

    return documented_fixture


def check_file_leaks(func) -> Callable:
    """
    Decorate a test function tot check that we are not leaking file descriptors.
    """
    psutil = safe_import("psutil")
    if not psutil:
        return func

    @wraps(func)
    def new_func(*args, **kwargs):
        proc = psutil.Process()
        flist = proc.open_files()

        func(*args, **kwargs)

        flist2 = proc.open_files()
        assert flist2 == flist

    return new_func


def async_mark():
    try:
        import_optional_dependency("pytest_asyncio")
        async_mark = pytest.mark.asyncio
    except ImportError:
        async_mark = pytest.mark.skip(reason="Missing dependency pytest-asyncio")

    return async_mark
