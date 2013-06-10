# pylint: disable-msg=W0614,W0401,W0611,W0622

__docformat__ = 'restructuredtext'

try:
    from pandas import hashtable, tslib, lib
except ImportError as e:  # pragma: no cover
    module = str(e).lstrip('cannot import name ')  # hack but overkill to use re
    raise ImportError("C extensions: {0} not built".format(module))

from datetime import datetime
import numpy as np

from pandas.version import version as __version__
from pandas.info import __doc__

# let init-time option registration happen
import pandas.core.config_init

from pandas.core.api import *
from pandas.sparse.api import *
from pandas.stats.api import *
from pandas.tseries.api import *
from pandas.io.api import *

from pandas.util.testing import debug

from pandas.tools.describe import value_range
from pandas.tools.merge import merge, concat, ordered_merge
from pandas.tools.pivot import pivot_table, crosstab
from pandas.tools.plotting import scatter_matrix, plot_params
from pandas.tools.tile import cut, qcut
from pandas.core.reshape import melt
