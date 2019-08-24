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
import locale
from typing import Optional

from _pytest.mark.structures import MarkDecorator
import pytest

from pandas.compat import is_platform_32bit, is_platform_windows
from pandas.compat.numpy import _np_version

from pandas.core.computation.expressions import _NUMEXPR_INSTALLED, _USE_NUMEXPR


def safe_import(mod_name, min_version=None):
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


def _skip_if_no_scipy():
    return not (
        safe_import("scipy.stats")
        and safe_import("scipy.sparse")
        and safe_import("scipy.interpolate")
        and safe_import("scipy.signal")
    )


def skip_if_installed(package: str,) -> MarkDecorator:
    """
    Skip a test if a package is installed.

    Parameters
    ----------
    package : str
        The name of the package.
    """
    return pytest.mark.skipif(
        safe_import(package), reason="Skipping because {} is installed.".format(package)
    )


def skip_if_no(package: str, min_version: Optional[str] = None) -> MarkDecorator:
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
    msg = "Could not import '{}'".format(package)
    if min_version:
        msg += " satisfying a min_version of {}".format(min_version)
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
    _skip_if_has_locale(),
    reason="Specific locale is set {lang}".format(lang=locale.getlocale()[0]),
)
skip_if_not_us_locale = pytest.mark.skipif(
    _skip_if_not_us_locale(),
    reason="Specific locale is set " "{lang}".format(lang=locale.getlocale()[0]),
)
skip_if_no_scipy = pytest.mark.skipif(
    _skip_if_no_scipy(), reason="Missing SciPy requirement"
)
skip_if_no_ne = pytest.mark.skipif(
    not _USE_NUMEXPR,
    reason="numexpr enabled->{enabled}, "
    "installed->{installed}".format(enabled=_USE_NUMEXPR, installed=_NUMEXPR_INSTALLED),
)


def skip_if_np_lt(ver_str, reason=None, *args, **kwds):
    if reason is None:
        reason = "NumPy %s or greater required" % ver_str
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

    Parameters:
    ----------
        args: iterable
            Positional arguments for docstring.

    Returns:
    -------
    documented_fixture: function
        The decorated function wrapped within a pytest
        ``parametrize_fixture_doc`` mark
    """

    def documented_fixture(fixture):
        fixture.__doc__ = fixture.__doc__.format(*args)
        return fixture

    return documented_fixture
