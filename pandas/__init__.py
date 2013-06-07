# pylint: disable-msg=W0614,W0401,W0611,W0622

__docformat__ = 'restructuredtext'

try:
    from . import hashtable, tslib, lib
except Exception:  # pragma: no cover
    import sys
    e = sys.exc_info()[1]  # Py25 and Py3 current exception syntax conflict
    print e
    if 'No module named lib' in str(e):
        raise ImportError('C extensions not built: if you installed already '
                          'verify that you are not importing from the source '
                          'directory')
    else:
        raise

from datetime import datetime
import numpy as np

from pandas.version import version as __version__
from pandas.info import __doc__

# let init-time option registration happen
import pandas.core.config_init

from pandas.core.api import (
        factorize, match, unique, value_counts, isnull, notnull, save, load,
        Categorical, Factor, set_printoptions, reset_printoptions,
        set_eng_float_format, Index, Int64Index, MultiIndex, Series, TimeSeries,
        DataFrame, Panel, Panel4D, groupby, pivot, get_dummies, lreshape, WidePanel,
        DateOffset, to_datetime, DatetimeIndex, Timestamp, date_range, bdate_range,
        Period, PeriodIndex, datetools, get_option, set_option, reset_option,
        describe_option, options, DateRange # deprecated
        )
from pandas.sparse.api import (
        SparseArray, SparseList, SparseSeries, SparseTimeSeries,
        SparseDataFrame, SparsePanel
        )
from pandas.stats.api import (
        ols, fama_macbeth, rolling_count, rolling_max, rolling_min,
        rolling_sum, rolling_mean, rolling_std, rolling_cov, rolling_corr,
        rolling_var, rolling_skew, rolling_kurt, rolling_quantile,
        rolling_median, rolling_apply, rolling_corr_pairwise, rolling_window,
        ewma, ewmvar, ewmstd, ewmvol, ewmcorr, ewmcov, expanding_count,
        expanding_max, expanding_min, expanding_sum, expanding_mean,
        expanding_std, expanding_cov, expanding_corr, expanding_var,
        expanding_skew, expanding_kurt, expanding_quantile, expanding_median,
        expanding_apply, expanding_corr_pairwise
        )
from pandas.tseries.api import (
    DatetimeIndex, date_range, bdate_range, infer_freq, Period, PeriodIndex,
    period_range, pnow, TimeGrouper, NaT, offsets
    )
from pandas.io.api import (
        read_csv, read_table, read_clipboard, read_fwf, to_clipboard, ExcelFile,
        ExcelWriter, read_excel, HDFStore, Term, get_store, read_hdf, read_html,
        read_sql, read_stata
        )

from pandas.util.testing import debug

from pandas.tools.describe import value_range
from pandas.tools.merge import merge, concat, ordered_merge
from pandas.tools.pivot import pivot_table, crosstab
from pandas.tools.plotting import scatter_matrix, plot_params
from pandas.tools.tile import cut, qcut
from pandas.core.reshape import melt

# import these so we can add them to __all__
import pandas.core.api as core_api
import pandas.sparse.api as sparse_api
import pandas.stats.api as stats_api
import pandas.tseries.api as tseries_api
import pandas.io.api as io_api
__all__ = ["debug", "value_range", "merge", "concat", "ordered_merge",
        "pivot_table", "crosstab", "scatter_matrix", "plot_params",
        "cut", "qcut", "melt"]
__all__.extend(core_api.__all__)
__all__.extend(sparse_api.__all__)
__all__.extend(stats_api.__all__)
__all__.extend(tseries_api.__all__)
__all__.extend(io_api.__all__)
