from pandas import *
import pandas as pd
from datetime import timedelta
from numpy.random import randn
from numpy.random import randint
from numpy.random import permutation
import pandas.util.testing as tm
import random
import numpy as np
try:
    from pandas.compat import range
except ImportError:
    pass

np.random.seed(1234)
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
