# flake8: noqa

import warnings
warnings.warn("The pandas.tslib module is deprecated and will be "
              "removed in a future version. Please import from "
              "the pandas or pandas.errors instead", FutureWarning, stacklevel=2)
from pandas._libs.tslib import (Timestamp, Timedelta,
                               NaT, OutOfBoundsDatetime)
