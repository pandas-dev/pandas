# pylint: disable-msg=W0614,W0401,W0611,W0622

__docformat__ = 'restructuredtext'

try:
    from . import hashtable, tslib, lib
except Exception:  # pragma: no cover
    import sys
    e = sys.exc_info()[1]  # Py25 and Py3 current exception syntax conflict
    print(e)
    if 'No module named lib' in str(e):
        raise ImportError('C extensions not built: if you installed already '
                          'verify that you are not importing from the source '
                          'directory')
    else:
        raise

from datetime import datetime
import numpy as np

# XXX: HACK for NumPy 1.5.1 to suppress warnings
try:
    np.seterr(all='ignore')
except Exception:  # pragma: no cover
    pass

# numpy versioning
from distutils.version import LooseVersion
_np_version = np.version.short_version
_np_version_under1p6 = LooseVersion(_np_version) < '1.6'
_np_version_under1p7 = LooseVersion(_np_version) < '1.7'
_np_version_under1p8 = LooseVersion(_np_version) < '1.8'

from pandas.version import version as __version__
from pandas.info import __doc__

# let init-time option registration happen
import pandas.core.config_init

from pandas.core.api import *
from pandas.sparse.api import *
from pandas.stats.api import *
from pandas.tseries.api import *
from pandas.io.api import *
from pandas.computation.api import *

from pandas.tools.describe import value_range
from pandas.tools.merge import merge, concat, ordered_merge
from pandas.tools.pivot import pivot_table, crosstab
from pandas.tools.plotting import scatter_matrix, plot_params
from pandas.tools.tile import cut, qcut
from pandas.core.reshape import melt
from pandas.util.print_versions import show_versions
