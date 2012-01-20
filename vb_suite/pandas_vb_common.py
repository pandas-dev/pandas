from pandas import *
from pandas.util.testing import rands
from datetime import timedelta
from numpy.random import randn
import pandas.util.testing as tm
import random
import numpy as np

try:
    Panel = WidePanel
except Exception:
    pass

# didn't add to namespace until later
try:
    from pandas.core.index import MultiIndex
except ImportError:
    pass
