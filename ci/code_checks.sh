#!/bin/bash
#
# Run checks related to code quality.
#
# This script is intended for both the CI and to check locally that code standards are
# respected. We run doctests here (currently some files only), and we
# validate formatting error in docstrings.
#
# Usage:
#   $ ./ci/code_checks.sh               # run all checks
#   $ ./ci/code_checks.sh code          # checks on imported code
#   $ ./ci/code_checks.sh doctests      # run doctests
#   $ ./ci/code_checks.sh docstrings    # validate docstring errors
#   $ ./ci/code_checks.sh single-docs   # check single-page docs build warning-free
#   $ ./ci/code_checks.sh notebooks     # check execution of documentation notebooks

set -uo pipefail

if [[ -v 1 ]]; then
    CHECK=$1
else
    # script will fail if it uses an unset variable (i.e. $1 is not provided)
    CHECK=""
fi

[[ -z "$CHECK" || "$CHECK" == "code" || "$CHECK" == "doctests" || "$CHECK" == "docstrings" || "$CHECK" == "single-docs" || "$CHECK" == "notebooks" ]] || \
    { echo "Unknown command $1. Usage: $0 [code|doctests|docstrings|single-docs|notebooks]"; exit 9999; }

BASE_DIR="$(dirname $0)/.."
RET=0

### CODE ###
if [[ -z "$CHECK" || "$CHECK" == "code" ]]; then

    MSG='Check import. No warnings, and blocklist some optional dependencies' ; echo $MSG
    python -W error -c "
import sys
import pandas

blocklist = {'bs4', 'gcsfs', 'html5lib', 'http', 'ipython', 'jinja2', 'hypothesis',
             'lxml', 'matplotlib', 'openpyxl', 'py', 'pytest', 's3fs', 'scipy',
             'tables', 'urllib.request', 'xlrd', 'xlsxwriter'}

# GH#28227 for some of these check for top-level modules, while others are
#  more specific (e.g. urllib.request)
import_mods = set(m.split('.')[0] for m in sys.modules) | set(sys.modules)
mods = blocklist & import_mods
if mods:
    sys.stderr.write('err: pandas should not import: {}\n'.format(', '.join(mods)))
    sys.exit(len(mods))
    "
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCTESTS ###
if [[ -z "$CHECK" || "$CHECK" == "doctests" ]]; then

    MSG='Python and Cython Doctests' ; echo $MSG
    python -c 'import pandas as pd; pd.test(run_doctests=True)'
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCSTRINGS ###
if [[ -z "$CHECK" || "$CHECK" == "docstrings" ]]; then

    MSG='Validate docstrings (EX01, EX03, EX04, GL01, GL02, GL03, GL04, GL05, GL06, GL07, GL09, GL10, PR03, PR04, PR05, PR06, PR08, PR09, PR10, RT01, RT02, RT04, RT05, SA02, SA03, SA04, SS01, SS02, SS03, SS04, SS05, SS06)' ; echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX01,EX03,EX04,GL01,GL02,GL03,GL04,GL05,GL06,GL07,GL09,GL10,PR03,PR04,PR05,PR06,PR08,PR09,PR10,RT01,RT02,RT04,RT05,SA02,SA03,SA04,SS01,SS02,SS03,SS04,SS05,SS06
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Partially validate docstrings (PR02)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=PR02 --ignore_functions \
        pandas.Series.dt.to_period\
        pandas.Series.dt.tz_localize\
        pandas.Series.dt.tz_convert\
        pandas.Series.dt.strftime\
        pandas.Series.dt.round\
        pandas.Series.dt.floor\
        pandas.Series.dt.ceil\
        pandas.Series.dt.month_name\
        pandas.Series.dt.day_name\
        pandas.Series.str.wrap\
        pandas.Series.cat.rename_categories\
        pandas.Series.cat.reorder_categories\
        pandas.Series.cat.add_categories\
        pandas.Series.cat.remove_categories\
        pandas.Series.cat.set_categories\
        pandas.Series.plot\
        pandas.DataFrame.plot\
        pandas.tseries.offsets.DateOffset\
        pandas.tseries.offsets.BusinessDay\
        pandas.tseries.offsets.BDay\
        pandas.tseries.offsets.BusinessHour\
        pandas.tseries.offsets.CustomBusinessDay\
        pandas.tseries.offsets.CDay\
        pandas.tseries.offsets.CustomBusinessHour\
        pandas.tseries.offsets.MonthEnd\
        pandas.tseries.offsets.MonthBegin\
        pandas.tseries.offsets.BusinessMonthEnd\
        pandas.tseries.offsets.BMonthEnd\
        pandas.tseries.offsets.BusinessMonthBegin\
        pandas.tseries.offsets.BMonthBegin\
        pandas.tseries.offsets.CustomBusinessMonthEnd\
        pandas.tseries.offsets.CBMonthEnd\
        pandas.tseries.offsets.CustomBusinessMonthBegin\
        pandas.tseries.offsets.CBMonthBegin\
        pandas.tseries.offsets.SemiMonthEnd\
        pandas.tseries.offsets.SemiMonthBegin\
        pandas.tseries.offsets.Week\
        pandas.tseries.offsets.WeekOfMonth\
        pandas.tseries.offsets.LastWeekOfMonth\
        pandas.tseries.offsets.BQuarterEnd\
        pandas.tseries.offsets.BQuarterBegin\
        pandas.tseries.offsets.QuarterEnd\
        pandas.tseries.offsets.QuarterBegin\
        pandas.tseries.offsets.BYearEnd\
        pandas.tseries.offsets.BYearBegin\
        pandas.tseries.offsets.YearEnd\
        pandas.tseries.offsets.YearBegin\
        pandas.tseries.offsets.FY5253\
        pandas.tseries.offsets.FY5253Quarter\
        pandas.tseries.offsets.Easter\
        pandas.tseries.offsets.Day\
        pandas.tseries.offsets.Hour\
        pandas.tseries.offsets.Minute\
        pandas.tseries.offsets.Second\
        pandas.tseries.offsets.Milli\
        pandas.tseries.offsets.Micro\
        pandas.tseries.offsets.Nano\
        pandas.Timestamp.max\
        pandas.Timestamp.min\
        pandas.Timestamp.resolution\
        pandas.Timedelta.max\
        pandas.Timedelta.min\
        pandas.Timedelta.resolution\
        pandas.Interval\
        pandas.Grouper\
        pandas.core.groupby.DataFrameGroupBy.nth\
        pandas.core.groupby.DataFrameGroupBy.rolling\
        pandas.core.groupby.SeriesGroupBy.nth\
        pandas.core.groupby.SeriesGroupBy.rolling\
        pandas.core.groupby.DataFrameGroupBy.plot\
        pandas.core.groupby.SeriesGroupBy.plot # There should be no backslash in the final line, please keep this comment in the last ignored function
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Partially validate docstrings (PR01)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=PR01 --ignore_functions \
        pandas.Categorical\
        pandas.Categorical.__array__\
        pandas.CategoricalIndex.equals\
        pandas.CategoricalIndex.map\
        pandas.DataFrame.at_time\
        pandas.DataFrame.backfill\
        pandas.DataFrame.get\
        pandas.DataFrame.pad\
        pandas.DataFrame.sem\
        pandas.DataFrame.sparse\
        pandas.DataFrame.std\
        pandas.DataFrame.swapaxes\
        pandas.DataFrame.var\
        pandas.DatetimeIndex.indexer_at_time\
        pandas.DatetimeIndex.snap\
        pandas.DatetimeIndex.std\
        pandas.ExcelFile\
        pandas.ExcelFile.parse\
        pandas.Grouper\
        pandas.HDFStore.append\
        pandas.HDFStore.put\
        pandas.Index.get_indexer_for\
        pandas.Index.identical\
        pandas.Index.putmask\
        pandas.Index.ravel\
        pandas.Index.str\
        pandas.Index.take\
        pandas.IntervalDtype\
        pandas.MultiIndex\
        pandas.Period.strftime\
        pandas.RangeIndex.from_range\
        pandas.Series.at_time\
        pandas.Series.backfill\
        pandas.Series.cat.add_categories\
        pandas.Series.cat.as_ordered\
        pandas.Series.cat.as_unordered\
        pandas.Series.cat.remove_categories\
        pandas.Series.cat.remove_unused_categories\
        pandas.Series.cat.rename_categories\
        pandas.Series.cat.reorder_categories\
        pandas.Series.cat.set_categories\
        pandas.Series.dt.ceil\
        pandas.Series.dt.day_name\
        pandas.Series.dt.floor\
        pandas.Series.dt.month_name\
        pandas.Series.dt.normalize\
        pandas.Series.dt.round\
        pandas.Series.dt.strftime\
        pandas.Series.dt.to_period\
        pandas.Series.dt.total_seconds\
        pandas.Series.dt.tz_convert\
        pandas.Series.dt.tz_localize\
        pandas.Series.get\
        pandas.Series.pad\
        pandas.Series.sem\
        pandas.Series.sparse\
        pandas.Series.std\
        pandas.Series.str\
        pandas.Series.str.wrap\
        pandas.Series.var\
        pandas.Timedelta.to_numpy\
        pandas.TimedeltaIndex\
        pandas.Timestamp.combine\
        pandas.Timestamp.fromtimestamp\
        pandas.Timestamp.strptime\
        pandas.Timestamp.to_numpy\
        pandas.Timestamp.to_period\
        pandas.Timestamp.to_pydatetime\
        pandas.Timestamp.utcfromtimestamp\
        pandas.api.extensions.ExtensionArray._pad_or_backfill\
        pandas.api.extensions.ExtensionArray.interpolate\
        pandas.api.indexers.BaseIndexer\
        pandas.api.indexers.FixedForwardWindowIndexer\
        pandas.api.indexers.VariableOffsetWindowIndexer\
        pandas.api.types.is_bool\
        pandas.api.types.is_complex\
        pandas.api.types.is_float\
        pandas.api.types.is_hashable\
        pandas.api.types.is_integer\
        pandas.core.groupby.DataFrameGroupBy.cummax\
        pandas.core.groupby.DataFrameGroupBy.cummin\
        pandas.core.groupby.DataFrameGroupBy.cumprod\
        pandas.core.groupby.DataFrameGroupBy.cumsum\
        pandas.core.groupby.DataFrameGroupBy.filter\
        pandas.core.groupby.DataFrameGroupBy.pct_change\
        pandas.core.groupby.DataFrameGroupBy.rolling\
        pandas.core.groupby.SeriesGroupBy.cummax\
        pandas.core.groupby.SeriesGroupBy.cummin\
        pandas.core.groupby.SeriesGroupBy.cumprod\
        pandas.core.groupby.SeriesGroupBy.cumsum\
        pandas.core.groupby.SeriesGroupBy.filter\
        pandas.core.groupby.SeriesGroupBy.nunique\
        pandas.core.groupby.SeriesGroupBy.pct_change\
        pandas.core.groupby.SeriesGroupBy.rolling\
        pandas.core.resample.Resampler.max\
        pandas.core.resample.Resampler.min\
        pandas.core.resample.Resampler.quantile\
        pandas.core.resample.Resampler.transform\
        pandas.core.window.expanding.Expanding.corr\
        pandas.core.window.expanding.Expanding.count\
        pandas.core.window.rolling.Rolling.max\
        pandas.core.window.rolling.Window.std\
        pandas.core.window.rolling.Window.var\
        pandas.errors.AbstractMethodError\
        pandas.errors.UndefinedVariableError\
        pandas.get_option\
        pandas.io.formats.style.Styler.to_excel\
        pandas.melt\
        pandas.option_context\
        pandas.read_fwf\
        pandas.reset_option\
        pandas.tseries.offsets.BQuarterBegin.is_month_end\
        pandas.tseries.offsets.BQuarterBegin.is_month_start\
        pandas.tseries.offsets.BQuarterBegin.is_quarter_end\
        pandas.tseries.offsets.BQuarterBegin.is_quarter_start\
        pandas.tseries.offsets.BQuarterBegin.is_year_end\
        pandas.tseries.offsets.BQuarterBegin.is_year_start\
        pandas.tseries.offsets.BQuarterEnd.is_month_end\
        pandas.tseries.offsets.BQuarterEnd.is_month_start\
        pandas.tseries.offsets.BQuarterEnd.is_quarter_end\
        pandas.tseries.offsets.BQuarterEnd.is_quarter_start\
        pandas.tseries.offsets.BQuarterEnd.is_year_end\
        pandas.tseries.offsets.BQuarterEnd.is_year_start\
        pandas.tseries.offsets.BYearBegin.is_month_end\
        pandas.tseries.offsets.BYearBegin.is_month_start\
        pandas.tseries.offsets.BYearBegin.is_quarter_end\
        pandas.tseries.offsets.BYearBegin.is_quarter_start\
        pandas.tseries.offsets.BYearBegin.is_year_end\
        pandas.tseries.offsets.BYearBegin.is_year_start\
        pandas.tseries.offsets.BYearEnd.is_month_end\
        pandas.tseries.offsets.BYearEnd.is_month_start\
        pandas.tseries.offsets.BYearEnd.is_quarter_end\
        pandas.tseries.offsets.BYearEnd.is_quarter_start\
        pandas.tseries.offsets.BYearEnd.is_year_end\
        pandas.tseries.offsets.BYearEnd.is_year_start\
        pandas.tseries.offsets.BusinessDay.is_month_end\
        pandas.tseries.offsets.BusinessDay.is_month_start\
        pandas.tseries.offsets.BusinessDay.is_quarter_end\
        pandas.tseries.offsets.BusinessDay.is_quarter_start\
        pandas.tseries.offsets.BusinessDay.is_year_end\
        pandas.tseries.offsets.BusinessDay.is_year_start\
        pandas.tseries.offsets.BusinessHour.is_month_end\
        pandas.tseries.offsets.BusinessHour.is_month_start\
        pandas.tseries.offsets.BusinessHour.is_quarter_end\
        pandas.tseries.offsets.BusinessHour.is_quarter_start\
        pandas.tseries.offsets.BusinessHour.is_year_end\
        pandas.tseries.offsets.BusinessHour.is_year_start\
        pandas.tseries.offsets.BusinessMonthBegin.is_month_end\
        pandas.tseries.offsets.BusinessMonthBegin.is_month_start\
        pandas.tseries.offsets.BusinessMonthBegin.is_quarter_end\
        pandas.tseries.offsets.BusinessMonthBegin.is_quarter_start\
        pandas.tseries.offsets.BusinessMonthBegin.is_year_end\
        pandas.tseries.offsets.BusinessMonthBegin.is_year_start\
        pandas.tseries.offsets.BusinessMonthEnd.is_month_end\
        pandas.tseries.offsets.BusinessMonthEnd.is_month_start\
        pandas.tseries.offsets.BusinessMonthEnd.is_quarter_end\
        pandas.tseries.offsets.BusinessMonthEnd.is_quarter_start\
        pandas.tseries.offsets.BusinessMonthEnd.is_year_end\
        pandas.tseries.offsets.BusinessMonthEnd.is_year_start\
        pandas.tseries.offsets.CustomBusinessDay.is_month_end\
        pandas.tseries.offsets.CustomBusinessDay.is_month_start\
        pandas.tseries.offsets.CustomBusinessDay.is_quarter_end\
        pandas.tseries.offsets.CustomBusinessDay.is_quarter_start\
        pandas.tseries.offsets.CustomBusinessDay.is_year_end\
        pandas.tseries.offsets.CustomBusinessDay.is_year_start\
        pandas.tseries.offsets.CustomBusinessHour.is_month_end\
        pandas.tseries.offsets.CustomBusinessHour.is_month_start\
        pandas.tseries.offsets.CustomBusinessHour.is_quarter_end\
        pandas.tseries.offsets.CustomBusinessHour.is_quarter_start\
        pandas.tseries.offsets.CustomBusinessHour.is_year_end\
        pandas.tseries.offsets.CustomBusinessHour.is_year_start\
        pandas.tseries.offsets.CustomBusinessMonthBegin.is_month_end\
        pandas.tseries.offsets.CustomBusinessMonthBegin.is_month_start\
        pandas.tseries.offsets.CustomBusinessMonthBegin.is_quarter_end\
        pandas.tseries.offsets.CustomBusinessMonthBegin.is_quarter_start\
        pandas.tseries.offsets.CustomBusinessMonthBegin.is_year_end\
        pandas.tseries.offsets.CustomBusinessMonthBegin.is_year_start\
        pandas.tseries.offsets.CustomBusinessMonthEnd.is_month_end\
        pandas.tseries.offsets.CustomBusinessMonthEnd.is_month_start\
        pandas.tseries.offsets.CustomBusinessMonthEnd.is_quarter_end\
        pandas.tseries.offsets.CustomBusinessMonthEnd.is_quarter_start\
        pandas.tseries.offsets.CustomBusinessMonthEnd.is_year_end\
        pandas.tseries.offsets.CustomBusinessMonthEnd.is_year_start\
        pandas.tseries.offsets.DateOffset.is_month_end\
        pandas.tseries.offsets.DateOffset.is_month_start\
        pandas.tseries.offsets.DateOffset.is_quarter_end\
        pandas.tseries.offsets.DateOffset.is_quarter_start\
        pandas.tseries.offsets.DateOffset.is_year_end\
        pandas.tseries.offsets.DateOffset.is_year_start\
        pandas.tseries.offsets.Day.is_month_end\
        pandas.tseries.offsets.Day.is_month_start\
        pandas.tseries.offsets.Day.is_quarter_end\
        pandas.tseries.offsets.Day.is_quarter_start\
        pandas.tseries.offsets.Day.is_year_end\
        pandas.tseries.offsets.Day.is_year_start\
        pandas.tseries.offsets.Easter.is_month_end\
        pandas.tseries.offsets.Easter.is_month_start\
        pandas.tseries.offsets.Easter.is_quarter_end\
        pandas.tseries.offsets.Easter.is_quarter_start\
        pandas.tseries.offsets.Easter.is_year_end\
        pandas.tseries.offsets.Easter.is_year_start\
        pandas.tseries.offsets.FY5253.is_month_end\
        pandas.tseries.offsets.FY5253.is_month_start\
        pandas.tseries.offsets.FY5253.is_quarter_end\
        pandas.tseries.offsets.FY5253.is_quarter_start\
        pandas.tseries.offsets.FY5253.is_year_end\
        pandas.tseries.offsets.FY5253.is_year_start\
        pandas.tseries.offsets.FY5253Quarter.is_month_end\
        pandas.tseries.offsets.FY5253Quarter.is_month_start\
        pandas.tseries.offsets.FY5253Quarter.is_quarter_end\
        pandas.tseries.offsets.FY5253Quarter.is_quarter_start\
        pandas.tseries.offsets.FY5253Quarter.is_year_end\
        pandas.tseries.offsets.FY5253Quarter.is_year_start\
        pandas.tseries.offsets.Hour.is_month_end\
        pandas.tseries.offsets.Hour.is_month_start\
        pandas.tseries.offsets.Hour.is_quarter_end\
        pandas.tseries.offsets.Hour.is_quarter_start\
        pandas.tseries.offsets.Hour.is_year_end\
        pandas.tseries.offsets.Hour.is_year_start\
        pandas.tseries.offsets.LastWeekOfMonth.is_month_end\
        pandas.tseries.offsets.LastWeekOfMonth.is_month_start\
        pandas.tseries.offsets.LastWeekOfMonth.is_quarter_end\
        pandas.tseries.offsets.LastWeekOfMonth.is_quarter_start\
        pandas.tseries.offsets.LastWeekOfMonth.is_year_end\
        pandas.tseries.offsets.LastWeekOfMonth.is_year_start\
        pandas.tseries.offsets.Micro.is_month_end\
        pandas.tseries.offsets.Micro.is_month_start\
        pandas.tseries.offsets.Micro.is_quarter_end\
        pandas.tseries.offsets.Micro.is_quarter_start\
        pandas.tseries.offsets.Micro.is_year_end\
        pandas.tseries.offsets.Micro.is_year_start\
        pandas.tseries.offsets.Milli.is_month_end\
        pandas.tseries.offsets.Milli.is_month_start\
        pandas.tseries.offsets.Milli.is_quarter_end\
        pandas.tseries.offsets.Milli.is_quarter_start\
        pandas.tseries.offsets.Milli.is_year_end\
        pandas.tseries.offsets.Milli.is_year_start\
        pandas.tseries.offsets.Minute.is_month_end\
        pandas.tseries.offsets.Minute.is_month_start\
        pandas.tseries.offsets.Minute.is_quarter_end\
        pandas.tseries.offsets.Minute.is_quarter_start\
        pandas.tseries.offsets.Minute.is_year_end\
        pandas.tseries.offsets.Minute.is_year_start\
        pandas.tseries.offsets.MonthBegin.is_month_end\
        pandas.tseries.offsets.MonthBegin.is_month_start\
        pandas.tseries.offsets.MonthBegin.is_quarter_end\
        pandas.tseries.offsets.MonthBegin.is_quarter_start\
        pandas.tseries.offsets.MonthBegin.is_year_end\
        pandas.tseries.offsets.MonthBegin.is_year_start\
        pandas.tseries.offsets.MonthEnd.is_month_end\
        pandas.tseries.offsets.MonthEnd.is_month_start\
        pandas.tseries.offsets.MonthEnd.is_quarter_end\
        pandas.tseries.offsets.MonthEnd.is_quarter_start\
        pandas.tseries.offsets.MonthEnd.is_year_end\
        pandas.tseries.offsets.MonthEnd.is_year_start\
        pandas.tseries.offsets.Nano.is_month_end\
        pandas.tseries.offsets.Nano.is_month_start\
        pandas.tseries.offsets.Nano.is_quarter_end\
        pandas.tseries.offsets.Nano.is_quarter_start\
        pandas.tseries.offsets.Nano.is_year_end\
        pandas.tseries.offsets.Nano.is_year_start\
        pandas.tseries.offsets.QuarterBegin.is_month_end\
        pandas.tseries.offsets.QuarterBegin.is_month_start\
        pandas.tseries.offsets.QuarterBegin.is_quarter_end\
        pandas.tseries.offsets.QuarterBegin.is_quarter_start\
        pandas.tseries.offsets.QuarterBegin.is_year_end\
        pandas.tseries.offsets.QuarterBegin.is_year_start\
        pandas.tseries.offsets.QuarterEnd.is_month_end\
        pandas.tseries.offsets.QuarterEnd.is_month_start\
        pandas.tseries.offsets.QuarterEnd.is_quarter_end\
        pandas.tseries.offsets.QuarterEnd.is_quarter_start\
        pandas.tseries.offsets.QuarterEnd.is_year_end\
        pandas.tseries.offsets.QuarterEnd.is_year_start\
        pandas.tseries.offsets.Second.is_month_end\
        pandas.tseries.offsets.Second.is_month_start\
        pandas.tseries.offsets.Second.is_quarter_end\
        pandas.tseries.offsets.Second.is_quarter_start\
        pandas.tseries.offsets.Second.is_year_end\
        pandas.tseries.offsets.Second.is_year_start\
        pandas.tseries.offsets.SemiMonthBegin.is_month_end\
        pandas.tseries.offsets.SemiMonthBegin.is_month_start\
        pandas.tseries.offsets.SemiMonthBegin.is_quarter_end\
        pandas.tseries.offsets.SemiMonthBegin.is_quarter_start\
        pandas.tseries.offsets.SemiMonthBegin.is_year_end\
        pandas.tseries.offsets.SemiMonthBegin.is_year_start\
        pandas.tseries.offsets.SemiMonthEnd.is_month_end\
        pandas.tseries.offsets.SemiMonthEnd.is_month_start\
        pandas.tseries.offsets.SemiMonthEnd.is_quarter_end\
        pandas.tseries.offsets.SemiMonthEnd.is_quarter_start\
        pandas.tseries.offsets.SemiMonthEnd.is_year_end\
        pandas.tseries.offsets.SemiMonthEnd.is_year_start\
        pandas.tseries.offsets.Tick.is_month_end\
        pandas.tseries.offsets.Tick.is_month_start\
        pandas.tseries.offsets.Tick.is_quarter_end\
        pandas.tseries.offsets.Tick.is_quarter_start\
        pandas.tseries.offsets.Tick.is_year_end\
        pandas.tseries.offsets.Tick.is_year_start\
        pandas.tseries.offsets.Week.is_month_end\
        pandas.tseries.offsets.Week.is_month_start\
        pandas.tseries.offsets.Week.is_quarter_end\
        pandas.tseries.offsets.Week.is_quarter_start\
        pandas.tseries.offsets.Week.is_year_end\
        pandas.tseries.offsets.Week.is_year_start\
        pandas.tseries.offsets.WeekOfMonth.is_month_end\
        pandas.tseries.offsets.WeekOfMonth.is_month_start\
        pandas.tseries.offsets.WeekOfMonth.is_quarter_end\
        pandas.tseries.offsets.WeekOfMonth.is_quarter_start\
        pandas.tseries.offsets.WeekOfMonth.is_year_end\
        pandas.tseries.offsets.WeekOfMonth.is_year_start\
        pandas.tseries.offsets.YearBegin.is_month_end\
        pandas.tseries.offsets.YearBegin.is_month_start\
        pandas.tseries.offsets.YearBegin.is_quarter_end\
        pandas.tseries.offsets.YearBegin.is_quarter_start\
        pandas.tseries.offsets.YearBegin.is_year_end\
        pandas.tseries.offsets.YearBegin.is_year_start\
        pandas.tseries.offsets.YearEnd.is_month_end\
        pandas.tseries.offsets.YearEnd.is_month_start\
        pandas.tseries.offsets.YearEnd.is_quarter_end\
        pandas.tseries.offsets.YearEnd.is_quarter_start\
        pandas.tseries.offsets.YearEnd.is_year_end\
        pandas.tseries.offsets.YearEnd.is_year_start # There should be no backslash in the final line, please keep this comment in the last ignored function
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCUMENTATION NOTEBOOKS ###
if [[ -z "$CHECK" || "$CHECK" == "notebooks" ]]; then

    MSG='Notebooks' ; echo $MSG
    jupyter nbconvert --execute $(find doc/source -name '*.ipynb') --to notebook
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### SINGLE-PAGE DOCS ###
if [[ -z "$CHECK" || "$CHECK" == "single-docs" ]]; then
    python doc/make.py --warnings-are-errors --no-browser --single pandas.Series.value_counts
    python doc/make.py --warnings-are-errors --no-browser --single pandas.Series.str.split
fi

exit $RET
