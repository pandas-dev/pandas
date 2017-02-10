# pylint: disable-msg=W0614,W0401,W0611,W0622

# flake8: noqa

__docformat__ = 'restructuredtext'

# Let users know if they're missing any of our hard dependencies
hard_dependencies = ("numpy", "pytz", "dateutil")
missing_dependencies = []

for dependency in hard_dependencies:
    try:
        __import__(dependency)
    except ImportError as e:
        missing_dependencies.append(dependency)

if missing_dependencies:
    raise ImportError(
        "Missing required dependencies {0}".format(missing_dependencies))
del hard_dependencies, dependency, missing_dependencies

# numpy compat
from pandas.compat.numpy import *

try:
    from pandas import hashtable, tslib, lib
except ImportError as e:  # pragma: no cover
    # hack but overkill to use re
    module = str(e).lstrip('cannot import name ')
    raise ImportError("C extension: {0} not built. If you want to import "
                      "pandas from the source directory, you may need to run "
                      "'python setup.py build_ext --inplace --force' to build "
                      "the C extensions first.".format(module))

from datetime import datetime
from pandas.info import __doc__

# let init-time option registration happen
import pandas.core.config_init

from pandas.core.api import *
from pandas.sparse.api import *
from pandas.stats.api import *
from pandas.tseries.api import *
from pandas.computation.api import *

from pandas.tools.concat import concat
from pandas.tools.merge import (merge, ordered_merge,
                                merge_ordered, merge_asof)
from pandas.tools.pivot import pivot_table, crosstab
from pandas.tools.plotting import scatter_matrix, plot_params
from pandas.tools.tile import cut, qcut
from pandas.tools.util import to_numeric
from pandas.core.reshape import melt
from pandas.util.print_versions import show_versions

from pandas.io.api import *

from pandas.util._tester import test

# use the closest tagged version if possible
from ._version import get_versions
v = get_versions()
__version__ = v.get('closest-tag', v['version'])
del get_versions, v
