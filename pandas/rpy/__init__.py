
# GH9602
# deprecate rpy to instead directly use rpy2

# flake8: noqa

import warnings
warnings.warn("The pandas.rpy module is deprecated and will be "
              "removed in a future version. We refer to external packages "
              "like rpy2. "
              "\nSee here for a guide on how to port your code to rpy2: "
              "http://pandas.pydata.org/pandas-docs/stable/r_interface.html",
              FutureWarning, stacklevel=2)

try:
    from .common import importr, r, load_data
except ImportError:
    pass
