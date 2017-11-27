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
class Foo():
    ...

Or individual functions:

@td.skip_if_some_reason
def test_foo():
    ...

For more information, refer to the ``pytest`` documentation on ``skipif``.
"""

import pytest


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
        version = getattr(sys.modules[mod_name], '__version__')
        if version:
            from distutils.version import LooseVersion
            if LooseVersion(version) >= LooseVersion(min_version):
                return mod

    return False


def _skip_if_no_mpl():
    mod = safe_import("matplotlib")
    if mod:
        mod.use("Agg", warn=False)
    else:
        return True


skip_if_no_mpl = pytest.mark.skipif(_skip_if_no_mpl(),
                                    reason="Missing matplotlib dependency")
