from pandas import *
from pandas.util.testing import rands
from datetime import timedelta
from numpy.random import randn
from numpy.random import randint
from numpy.random import permutation
import pandas.util.compat as compat
import pandas.util.testing as tm
import random
import numpy as np

try:
    import pandas._tseries as lib
except:
    import pandas.lib as lib

try:
    Panel = WidePanel
except Exception:
    pass

# didn't add to namespace until later
try:
    from pandas.core.index import MultiIndex
except ImportError:
    pass
try:
    # if no range in py3compat, then don't import zip or map either
    from pandas.util.py3compat import range
    from six.moves import zip, map
except ImportError:
    pass
