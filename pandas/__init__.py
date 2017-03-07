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
    from pandas._libs import (hashtable as _hashtable,
                             lib as _lib,
                             tslib as _tslib)
except ImportError as e:  # pragma: no cover
    # hack but overkill to use re
    module = str(e).lstrip('cannot import name ')
    raise ImportError("C extension: {0} not built. If you want to import "
                      "pandas from the source directory, you may need to run "
                      "'python setup.py build_ext --inplace --force' to build "
                      "the C extensions first.".format(module))

from datetime import datetime

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

# extension module deprecations
from pandas.util.depr_module import _DeprecatedModule

json = _DeprecatedModule(deprmod='pandas.json', deprmodto='pandas.io.json.libjson')
parser = _DeprecatedModule(deprmod='pandas.parser', deprmodto='pandas.io.libparsers')
lib = _DeprecatedModule(deprmod='pandas.lib', deprmodto='pandas._libs.lib')
tslib = _DeprecatedModule(deprmod='pandas.tslib', deprmodto='pandas._libs.tslib')

# use the closest tagged version if possible
from ._version import get_versions
v = get_versions()
__version__ = v.get('closest-tag', v['version'])
del get_versions, v

# module level doc-string
__doc__ = """
pandas - a powerful data analysis and manipulation library for Python
=====================================================================

**pandas** is a Python package providing fast, flexible, and expressive data
structures designed to make working with "relational" or "labeled" data both
easy and intuitive. It aims to be the fundamental high-level building block for
doing practical, **real world** data analysis in Python. Additionally, it has
the broader goal of becoming **the most powerful and flexible open source data
analysis / manipulation tool available in any language**. It is already well on
its way toward this goal.

Main Features
-------------
Here are just a few of the things that pandas does well:

  - Easy handling of missing data in floating point as well as non-floating
    point data
  - Size mutability: columns can be inserted and deleted from DataFrame and
    higher dimensional objects
  - Automatic and explicit data alignment: objects can  be explicitly aligned
    to a set of labels, or the user can simply ignore the labels and let
    `Series`, `DataFrame`, etc. automatically align the data for you in
    computations
  - Powerful, flexible group by functionality to perform split-apply-combine
    operations on data sets, for both aggregating and transforming data
  - Make it easy to convert ragged, differently-indexed data in other Python
    and NumPy data structures into DataFrame objects
  - Intelligent label-based slicing, fancy indexing, and subsetting of large
    data sets
  - Intuitive merging and joining data sets
  - Flexible reshaping and pivoting of data sets
  - Hierarchical labeling of axes (possible to have multiple labels per tick)
  - Robust IO tools for loading data from flat files (CSV and delimited),
    Excel files, databases, and saving/loading data from the ultrafast HDF5
    format
  - Time series-specific functionality: date range generation and frequency
    conversion, moving window statistics, moving window linear regressions,
    date shifting and lagging, etc.
"""
