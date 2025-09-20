import os
import shutil
import sys
import tempfile
from importlib import import_module

import pytest

# Expose the unittest-driven decorators
from .ipunittest import ipdoctest, ipdocstring

def skipif(skip_condition, msg=None):
    """Make function raise SkipTest exception if skip_condition is true

    Parameters
    ----------

    skip_condition : bool or callable
      Flag to determine whether to skip test. If the condition is a
      callable, it is used at runtime to dynamically make the decision. This
      is useful for tests that may require costly imports, to delay the cost
      until the test suite is actually executed.
    msg : string
      Message to give on raising a SkipTest exception.

    Returns
    -------
    decorator : function
      Decorator, which, when applied to a function, causes SkipTest
      to be raised when the skip_condition was True, and the function
      to be called normally otherwise.
    """
    if msg is None:
        msg = "Test skipped due to test condition."

    assert isinstance(skip_condition, bool)
    return pytest.mark.skipif(skip_condition, reason=msg)


# A version with the condition set to true, common case just to attach a message
# to a skip decorator
def skip(msg=None):
    """Decorator factory - mark a test function for skipping from test suite.

    Parameters
    ----------
      msg : string
        Optional message to be added.

    Returns
    -------
       decorator : function
         Decorator, which, when applied to a function, causes SkipTest
         to be raised, with the optional message added.
    """
    if msg and not isinstance(msg, str):
        raise ValueError(
            "invalid object passed to `@skip` decorator, did you "
            "meant `@skip()` with brackets ?"
        )
    return skipif(True, msg)


def onlyif(condition, msg):
    """The reverse from skipif, see skipif for details."""

    return skipif(not condition, msg)


# -----------------------------------------------------------------------------
# Utility functions for decorators
def module_not_available(module):
    """Can module be imported?  Returns true if module does NOT import.

    This is used to make a decorator to skip tests that require module to be
    available, but delay the 'import numpy' to test execution time.
    """
    try:
        mod = import_module(module)
        mod_not_avail = False
    except ImportError:
        mod_not_avail = True

    return mod_not_avail


# -----------------------------------------------------------------------------
# Decorators for public use

# Decorators to skip certain tests on specific platforms.
skip_win32 = skipif(sys.platform == "win32", "This test does not run under Windows")


# Decorators to skip tests if not on specific platforms.
skip_if_not_win32 = skipif(sys.platform != "win32", "This test only runs under Windows")
skip_if_not_osx = skipif(
    not sys.platform.startswith("darwin"), "This test only runs under macOS"
)

_x11_skip_cond = (
    sys.platform not in ("darwin", "win32") and os.environ.get("DISPLAY", "") == ""
)
_x11_skip_msg = "Skipped under *nix when X11/XOrg not available"

skip_if_no_x11 = skipif(_x11_skip_cond, _x11_skip_msg)

# Other skip decorators

# generic skip without module
skip_without = lambda mod: skipif(
    module_not_available(mod), "This test requires %s" % mod
)

skipif_not_numpy = skip_without("numpy")

skipif_not_matplotlib = skip_without("matplotlib")

# A null 'decorator', useful to make more readable code that needs to pick
# between different decorators based on OS or other conditions
null_deco = lambda f: f

# Some tests only run where we can use unicode paths. Note that we can't just
# check os.path.supports_unicode_filenames, which is always False on Linux.
try:
    f = tempfile.NamedTemporaryFile(prefix="tmp€")
except UnicodeEncodeError:
    unicode_paths = False
# TODO: should this be finnally ?
else:
    unicode_paths = True
    f.close()

onlyif_unicode_paths = onlyif(
    unicode_paths,
    ("This test is only applicable where we can use unicode in filenames."),
)


def onlyif_cmds_exist(*commands):
    """
    Decorator to skip test when at least one of `commands` is not found.
    """
    for cmd in commands:
        reason = f"This test runs only if command '{cmd}' is installed"
        if not shutil.which(cmd):

            return pytest.mark.skip(reason=reason)
    return null_deco
