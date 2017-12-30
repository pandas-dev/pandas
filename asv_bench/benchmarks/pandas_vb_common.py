import os
from pandas import *
import pandas as pd
from numpy.random import randn
from numpy.random import randint
import pandas.util.testing as tm
import random
import numpy as np
import threading
from importlib import import_module

try:
    from pandas.compat import range
except ImportError:
    pass

numeric_dtypes = [np.int64, np.int32, np.uint32, np.uint64, np.float32,
                  np.float64, np.int16, np.int8, np.uint16, np.uint8]
datetime_dtypes = [np.datetime64, np.timedelta64]

# This function just needs to be imported into each benchmark file in order to
# sets up the random seed before each function.
# http://asv.readthedocs.io/en/latest/writing_benchmarks.html
def setup(*args, **kwargs):
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
        except:
            # On Windows, attempting to remove a file that is in use
            # causes an exception to be raised
            pass

    def teardown(self, *args, **kwargs):
        self.remove(self.fname)

# Compatibility import for lib
for imp in ['pandas._libs.lib', 'pandas.lib', 'pandas_tseries']:
    try:
        lib = import_module(imp)
        break
    except:
        pass

try:
    Panel = Panel
except Exception:
    Panel = WidePanel

# didn't add to namespace until later
try:
    from pandas.core.index import MultiIndex
except ImportError:
    pass
