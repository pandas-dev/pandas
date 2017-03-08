# flake8: noqa

import warnings
warnings.warn("The pandas.json module is deprecated and will be "
              "removed in a future version. Please import from "
              "the pandas.io.json instead", FutureWarning, stacklevel=2)
from pandas.io.json.libjson import dumps, loads
