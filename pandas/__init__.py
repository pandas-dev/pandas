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
    module = str(e).replace('cannot import name ', '')
    raise ImportError("C extension: {0} not built. If you want to import "
                      "pandas from the source directory, you may need to run "
                      "'python setup.py build_ext --inplace --force' to build "
                      "the C extensions first.".format(module))

from datetime import datetime

from pandas._config import (get_option, set_option, reset_option,
                            describe_option, option_context, options)

# let init-time option registration happen
import pandas.core.config_init

from pandas.core.panel import Panel
from pandas.core.arrays.sparse import SparseArray, SparseDtype
from pandas.core.sparse.frame import SparseDataFrame
from pandas.core.sparse.series import SparseSeries
from pandas.core.computation.eval import eval
from pandas.core.reshape.concat import concat
from pandas.core.reshape.melt import lreshape, melt, wide_to_long
from pandas.core.reshape.merge import merge, merge_asof, merge_ordered
from pandas.core.reshape.pivot import crosstab, pivot, pivot_table
from pandas.core.reshape.reshape import get_dummies
from pandas.core.reshape.tile import cut, qcut
from pandas.core.internals import BlockManager, _safe_reshape, make_block

from pandas.core.index import (
    Float64Index, Index, DatetimeIndex, Int64Index, MultiIndex, PeriodIndex, RangeIndex,
    TimedeltaIndex, CategoricalIndex)
from pandas.core.arrays import DatetimeArray, IntervalArray, PeriodArray
from pandas.core.arrays.sparse import BlockIndex, IntIndex
from pandas.core.arrays.categorical import Categorical
from pandas.core.arrays.interval import Interval
from pandas.core.indexes.interval import IntervalIndex
from pandas.core.indexes.datetimes import DatetimeIndex
from pandas.core.missing import isna
from pandas.core.tools.datetimes import to_datetime
from pandas.core.tools.timedeltas import to_timedelta
from pandas._libs.tslibs.timedeltas import Timedelta
from pandas._libs.tslibs.timestamps import Timestamp
from pandas._libs import NaT, Period, join
from pandas.io.clipboards import read_clipboard
from pandas.io.excel import ExcelFile, ExcelWriter, read_excel
from pandas.io.feather_format import read_feather
from pandas.io.gbq import read_gbq
from pandas.io.html import read_html
from pandas.io.json import read_json
from pandas.io.packers import read_msgpack, to_msgpack
from pandas.io.parquet import read_parquet
from pandas.io.parsers import read_csv, read_fwf, read_table
from pandas.io.pickle import read_pickle, to_pickle
from pandas.io.pytables import HDFStore, read_hdf
from pandas.io.sas import read_sas
from pandas.io.sql import read_sql, read_sql_query, read_sql_table
from pandas.io.stata import read_stata
from pandas.io.msgpack import ExtType, Packer as _Packer, Unpacker as _Unpacker
from pandas.io.formats.printing import pprint_thing
from pandas.io.parsers import TextParser
from pandas.io.common import _stringify_path, get_filepath_or_buffer, is_s3_url

from pandas.tseries.frequencies import infer_freq
from pandas.tseries import offsets

from pandas.util._print_versions import show_versions
from pandas.util._tester import test
import pandas.testing
import pandas.arrays

# use the closest tagged version if possible
from ._version import get_versions
v = get_versions()
__version__ = v.get('closest-tag', v['version'])
__git_version__ = v.get('full-revisionid')
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
    point data.
  - Size mutability: columns can be inserted and deleted from DataFrame and
    higher dimensional objects
  - Automatic and explicit data alignment: objects can be explicitly aligned
    to a set of labels, or the user can simply ignore the labels and let
    `Series`, `DataFrame`, etc. automatically align the data for you in
    computations.
  - Powerful, flexible group by functionality to perform split-apply-combine
    operations on data sets, for both aggregating and transforming data.
  - Make it easy to convert ragged, differently-indexed data in other Python
    and NumPy data structures into DataFrame objects.
  - Intelligent label-based slicing, fancy indexing, and subsetting of large
    data sets.
  - Intuitive merging and joining data sets.
  - Flexible reshaping and pivoting of data sets.
  - Hierarchical labeling of axes (possible to have multiple labels per tick).
  - Robust IO tools for loading data from flat files (CSV and delimited),
    Excel files, databases, and saving/loading data from the ultrafast HDF5
    format.
  - Time series-specific functionality: date range generation and frequency
    conversion, moving window statistics, moving window linear regressions,
    date shifting and lagging, etc.
"""
