import os
from importlib import import_module

import numpy as np
from pandas import Panel

# Compatibility import for lib
for imp in ['pandas._libs.lib', 'pandas.lib']:
    try:
        lib = import_module(imp)
        break
    except (ImportError, TypeError, ValueError):
        pass

numeric_dtypes = [np.int64, np.int32, np.uint32, np.uint64, np.float32,
                  np.float64, np.int16, np.int8, np.uint16, np.uint8]
datetime_dtypes = [np.datetime64, np.timedelta64]


def setup(*args, **kwargs):
    # This function just needs to be imported into each benchmark file to
    # set up the random seed before each function.
    # http://asv.readthedocs.io/en/latest/writing_benchmarks.html
    np.random.seed(1234)


class BaseIO(object):
    """
    Base class for IO benchmarks
    """
    fname = None

    def remove(self, f):
        """Remove created files"""
        try:
            os.remove(f)
        except OSError:
            # On Windows, attempting to remove a file that is in use
            # causes an exception to be raised
            pass

    def teardown(self, *args, **kwargs):
        self.remove(self.fname)
