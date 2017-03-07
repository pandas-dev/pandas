# flake8: noqa

import warnings
warnings.warn("The pandas.lib module is deprecated and will be "
              "removed in a future version. Please import from "
              "the pandas._libs.lib instead", FutureWarning, stacklevel=2)
from pandas._libs.lib import *
