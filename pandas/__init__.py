# pylint: disable-msg=W0614,W0401,W0611,W0622


__docformat__ = 'restructuredtext'

try:
    from pandas import hashtable, tslib, lib
except ImportError as e:  # pragma: no cover
    module = str(e).lstrip('cannot import name ')  # hack but overkill to use re
    raise ImportError("C extension: {0} not built. If you want to import "
                      "pandas from the source directory, you may need to run "
                      "'python setup.py build_ext --inplace' to build the C "
                      "extensions first.".format(module))

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
_np_version_under1p8 = LooseVersion(_np_version) < '1.8'
_np_version_under1p9 = LooseVersion(_np_version) < '1.9'


from pandas.info import __doc__


if LooseVersion(_np_version) < '1.7.0':
    raise ImportError('pandas {0} is incompatible with numpy < 1.7.0, '
                      'your numpy version is {1}. Please upgrade numpy to'
                      ' >= 1.7.0 to use pandas version {0}'.format(__version__,
                                                                   _np_version))

# let init-time option registration happen
import pandas.core.config_init

from pandas.core.api import *
from pandas.sparse.api import *
from pandas.stats.api import *
from pandas.tseries.api import *
from pandas.io.api import *
from pandas.computation.api import *

from pandas.tools.merge import merge, concat, ordered_merge
from pandas.tools.pivot import pivot_table, crosstab
from pandas.tools.plotting import scatter_matrix, plot_params
from pandas.tools.tile import cut, qcut
from pandas.tools.util import to_numeric
from pandas.core.reshape import melt
from pandas.util.print_versions import show_versions
import pandas.util.testing

# use the closest tagged version if possible
from ._version import get_versions
v = get_versions()
__version__ = v.get('closest-tag',v['version'])
del get_versions, v
