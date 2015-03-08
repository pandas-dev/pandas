
# GH9602
# deprecate rpy to instead directly use rpy2

import warnings
warnings.warn("The pandas.rpy module is deprecated and will be "
              "removed in a future version. We refer to external packages "
              "like rpy2, found here: http://rpy.sourceforge.net", FutureWarning)

try:
    from .common import importr, r, load_data
except ImportError:
    pass
