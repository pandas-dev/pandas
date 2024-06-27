# Licensed under a 3-clause BSD style license - see LICENSE.rst

import os
import sys

if sys.version_info[0] == 3:
    xrange = range


# Test that local imports work
from .shared import shared_function

# Test that asv's internal modules aren't visible on top level
if sys.version_info[0] < 3:
    import commands
try:
    import commands.quickstart  # noqa F401 unused import
    assert False
except ImportError:
    # OK
    pass


class TimeSecondary:
    sample_time = 0.05
    _printed = False

    def time_factorial(self):
        x = 1
        for i in xrange(100):
            x *= i
        # This is to print things to stdout, but not spam too much
        if not self._printed:
            sys.stdout.write("X")
            self._printed = True

    def time_exception(self):
        raise RuntimeError()


def track_value():
    return 42.0


def test_shared_code():
    assert shared_function() == 42


def track_environment_value():
    v = os.environ.get('SOME_TEST_VAR', '0')
    try:
        return int(v)
    except (ValueError, TypeError):
        return 0


def track_fail_errcode_123():
    if hasattr(os, '_exit'):
        os._exit(123)
    else:
        sys.exit(123)


def track_fail_signal_9():
    os.kill(os.getpid(), 9)
