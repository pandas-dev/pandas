from __future__ import annotations
__docformat__ = 'restructuredtext'
_hard_dependencies = ('numpy', 'dateutil')
for _dependency in _hard_dependencies:
    try:
        __import__(_dependency)
    except ImportError as _e:
        raise ImportError(f'Unable to import required dependency {_dependency}. Please see the traceback for details.') from _e
del _hard_dependencies, _dependency
try:
    from pandas.compat import is_numpy_dev as _is_numpy_dev
except ImportError as _err:
    _module = _err.name
    raise ImportError(f"C extension: {_module} not built. If you want to import pandas from the source directory, you may need to run 'python -m pip install -ve . --no-build-isolation -Ceditable-verbose=true' to build the C extensions first.") from _err
from pandas._config import get_option, set_option, reset_option, describe_option, option_context, options
import pandas.core.config_init
from pandas.core.api import ArrowDtype, Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype, UInt8Dtype, UInt16Dtype, UInt32Dtype, UInt64Dtype, Float32Dtype, Float64Dtype, CategoricalDtype, PeriodDtype, IntervalDtype, DatetimeTZDtype, StringDtype, BooleanDtype, NA, isna, isnull, notna, notnull, Index, CategoricalIndex, RangeIndex, MultiIndex, IntervalIndex, TimedeltaIndex, DatetimeIndex, PeriodIndex, IndexSlice, NaT, Period, period_range, Timedelta, timedelta_range, Timestamp, date_range, bdate_range, Interval, interval_range, DateOffset, to_numeric, to_datetime, to_timedelta, Flags, Grouper, factorize, unique, NamedAgg, NamedFunc, array, Categorical, set_eng_float_format, Series, DataFrame
from pandas.core.col import col
from pandas.core.dtypes.dtypes import SparseDtype
from pandas.tseries.api import infer_freq
from pandas.tseries import offsets
from pandas.core.computation.api import eval
from pandas.core.reshape.api import concat, lreshape, melt, wide_to_long, merge, merge_asof, merge_ordered, crosstab, pivot, pivot_table, get_dummies, from_dummies, cut, qcut
from pandas import api, arrays, errors, io, plotting, tseries
from pandas import testing
from pandas.util._print_versions import show_versions
from pandas.io.api import ExcelFile, ExcelWriter, read_excel, read_csv, read_fwf, read_table, read_pickle, to_pickle, HDFStore, read_hdf, read_sql, read_sql_query, read_sql_table, read_clipboard, read_parquet, read_orc, read_feather, read_html, read_xml, read_json, read_stata, read_sas, read_spss, read_iceberg
from pandas.io.json._normalize import json_normalize
from pandas.util._tester import test
from pandas._version_meson import __version__, __git_version__
__doc__ = '\npandas - a powerful data analysis and manipulation library for Python\n=====================================================================\n\n**pandas** is a Python package providing fast, flexible, and expressive data\nstructures designed to make working with "relational" or "labeled" data both\neasy and intuitive. It aims to be the fundamental high-level building block for\ndoing practical, **real world** data analysis in Python. Additionally, it has\nthe broader goal of becoming **the most powerful and flexible open source data\nanalysis / manipulation tool available in any language**. It is already well on\nits way toward this goal.\n\nMain Features\n-------------\nHere are just a few of the things that pandas does well:\n\n  - Easy handling of missing data in floating point as well as non-floating\n    point data.\n  - Size mutability: columns can be inserted and deleted from DataFrame and\n    higher dimensional objects\n  - Automatic and explicit data alignment: objects can be explicitly aligned\n    to a set of labels, or the user can simply ignore the labels and let\n    `Series`, `DataFrame`, etc. automatically align the data for you in\n    computations.\n  - Powerful, flexible group by functionality to perform split-apply-combine\n    operations on data sets, for both aggregating and transforming data.\n  - Make it easy to convert ragged, differently-indexed data in other Python\n    and NumPy data structures into DataFrame objects.\n  - Intelligent label-based slicing, fancy indexing, and subsetting of large\n    data sets.\n  - Intuitive merging and joining data sets.\n  - Flexible reshaping and pivoting of data sets.\n  - Hierarchical labeling of axes (possible to have multiple labels per tick).\n  - Robust IO tools for loading data from flat files (CSV and delimited),\n    Excel files, databases, and saving/loading data from the ultrafast HDF5\n    format.\n  - Time series-specific functionality: date range generation and frequency\n    conversion, moving window statistics, date shifting and lagging.\n'
__all__ = ['NA', 'ArrowDtype', 'BooleanDtype', 'Categorical', 'CategoricalDtype', 'CategoricalIndex', 'DataFrame', 'DateOffset', 'DatetimeIndex', 'DatetimeTZDtype', 'ExcelFile', 'ExcelWriter', 'Flags', 'Float32Dtype', 'Float64Dtype', 'Grouper', 'HDFStore', 'Index', 'IndexSlice', 'Int8Dtype', 'Int16Dtype', 'Int32Dtype', 'Int64Dtype', 'Interval', 'IntervalDtype', 'IntervalIndex', 'MultiIndex', 'NaT', 'NamedAgg', 'NamedFunc', 'Period', 'PeriodDtype', 'PeriodIndex', 'RangeIndex', 'Series', 'SparseDtype', 'StringDtype', 'Timedelta', 'TimedeltaIndex', 'Timestamp', 'UInt8Dtype', 'UInt16Dtype', 'UInt32Dtype', 'UInt64Dtype', '__git_version__', '__version__', 'api', 'array', 'arrays', 'bdate_range', 'col', 'concat', 'crosstab', 'cut', 'date_range', 'describe_option', 'errors', 'eval', 'factorize', 'from_dummies', 'get_dummies', 'get_option', 'infer_freq', 'interval_range', 'io', 'isna', 'isnull', 'json_normalize', 'lreshape', 'melt', 'merge', 'merge_asof', 'merge_ordered', 'notna', 'notnull', 'offsets', 'option_context', 'options', 'period_range', 'pivot', 'pivot_table', 'plotting', 'qcut', 'read_clipboard', 'read_csv', 'read_excel', 'read_feather', 'read_fwf', 'read_hdf', 'read_html', 'read_iceberg', 'read_json', 'read_orc', 'read_parquet', 'read_pickle', 'read_sas', 'read_spss', 'read_sql', 'read_sql_query', 'read_sql_table', 'read_stata', 'read_table', 'read_xml', 'reset_option', 'set_eng_float_format', 'set_option', 'show_versions', 'test', 'testing', 'timedelta_range', 'to_datetime', 'to_numeric', 'to_pickle', 'to_timedelta', 'tseries', 'unique', 'wide_to_long']
