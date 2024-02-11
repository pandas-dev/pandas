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

    MSG='Partially validate docstrings (GL08)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=GL08 --ignore_functions \
        pandas.DatetimeIndex.as_unit\
        pandas.DatetimeIndex.freq\
        pandas.ExcelFile.book\
        pandas.ExcelFile.sheet_names\
        pandas.Index.empty\
        pandas.Index.names\
        pandas.Index.view\
        pandas.IntervalIndex.left\
        pandas.IntervalIndex.length\
        pandas.IntervalIndex.mid\
        pandas.IntervalIndex.right\
        pandas.MultiIndex.codes\
        pandas.Period.freq\
        pandas.Period.ordinal\
        pandas.PeriodIndex.freq\
        pandas.PeriodIndex.qyear\
        pandas.Series.dt\
        pandas.Series.dt.as_unit\
        pandas.Series.dt.freq\
        pandas.Series.dt.qyear\
        pandas.Series.dt.unit\
        pandas.Series.empty\
        pandas.Timedelta.microseconds\
        pandas.Timedelta.unit\
        pandas.Timedelta.value\
        pandas.Timestamp.day\
        pandas.Timestamp.fold\
        pandas.Timestamp.hour\
        pandas.Timestamp.microsecond\
        pandas.Timestamp.minute\
        pandas.Timestamp.month\
        pandas.Timestamp.nanosecond\
        pandas.Timestamp.second\
        pandas.Timestamp.tzinfo\
        pandas.Timestamp.value\
        pandas.Timestamp.year\
        pandas.core.groupby.SeriesGroupBy.value_counts\
        pandas.tseries.offsets.BQuarterBegin.is_anchored\
        pandas.tseries.offsets.BQuarterBegin.is_on_offset\
        pandas.tseries.offsets.BQuarterBegin.n\
        pandas.tseries.offsets.BQuarterBegin.nanos\
        pandas.tseries.offsets.BQuarterBegin.normalize\
        pandas.tseries.offsets.BQuarterBegin.rule_code\
        pandas.tseries.offsets.BQuarterBegin.startingMonth\
        pandas.tseries.offsets.BQuarterEnd.is_anchored\
        pandas.tseries.offsets.BQuarterEnd.is_on_offset\
        pandas.tseries.offsets.BQuarterEnd.n\
        pandas.tseries.offsets.BQuarterEnd.nanos\
        pandas.tseries.offsets.BQuarterEnd.normalize\
        pandas.tseries.offsets.BQuarterEnd.rule_code\
        pandas.tseries.offsets.BQuarterEnd.startingMonth\
        pandas.tseries.offsets.BYearBegin.is_on_offset\
        pandas.tseries.offsets.BYearBegin.month\
        pandas.tseries.offsets.BYearBegin.n\
        pandas.tseries.offsets.BYearBegin.nanos\
        pandas.tseries.offsets.BYearBegin.normalize\
        pandas.tseries.offsets.BYearBegin.rule_code\
        pandas.tseries.offsets.BYearEnd.is_on_offset\
        pandas.tseries.offsets.BYearEnd.month\
        pandas.tseries.offsets.BYearEnd.n\
        pandas.tseries.offsets.BYearEnd.nanos\
        pandas.tseries.offsets.BYearEnd.normalize\
        pandas.tseries.offsets.BYearEnd.rule_code\
        pandas.tseries.offsets.BusinessDay.calendar\
        pandas.tseries.offsets.BusinessDay.holidays\
        pandas.tseries.offsets.BusinessDay.is_on_offset\
        pandas.tseries.offsets.BusinessDay.n\
        pandas.tseries.offsets.BusinessDay.nanos\
        pandas.tseries.offsets.BusinessDay.normalize\
        pandas.tseries.offsets.BusinessDay.rule_code\
        pandas.tseries.offsets.BusinessDay.weekmask\
        pandas.tseries.offsets.BusinessHour.calendar\
        pandas.tseries.offsets.BusinessHour.end\
        pandas.tseries.offsets.BusinessHour.holidays\
        pandas.tseries.offsets.BusinessHour.is_on_offset\
        pandas.tseries.offsets.BusinessHour.n\
        pandas.tseries.offsets.BusinessHour.nanos\
        pandas.tseries.offsets.BusinessHour.normalize\
        pandas.tseries.offsets.BusinessHour.rule_code\
        pandas.tseries.offsets.BusinessHour.start\
        pandas.tseries.offsets.BusinessHour.weekmask\
        pandas.tseries.offsets.BusinessMonthBegin.is_on_offset\
        pandas.tseries.offsets.BusinessMonthBegin.n\
        pandas.tseries.offsets.BusinessMonthBegin.nanos\
        pandas.tseries.offsets.BusinessMonthBegin.normalize\
        pandas.tseries.offsets.BusinessMonthBegin.rule_code\
        pandas.tseries.offsets.BusinessMonthEnd.is_on_offset\
        pandas.tseries.offsets.BusinessMonthEnd.n\
        pandas.tseries.offsets.BusinessMonthEnd.nanos\
        pandas.tseries.offsets.BusinessMonthEnd.normalize\
        pandas.tseries.offsets.BusinessMonthEnd.rule_code\
        pandas.tseries.offsets.CustomBusinessDay.calendar\
        pandas.tseries.offsets.CustomBusinessDay.holidays\
        pandas.tseries.offsets.CustomBusinessDay.is_on_offset\
        pandas.tseries.offsets.CustomBusinessDay.n\
        pandas.tseries.offsets.CustomBusinessDay.nanos\
        pandas.tseries.offsets.CustomBusinessDay.normalize\
        pandas.tseries.offsets.CustomBusinessDay.rule_code\
        pandas.tseries.offsets.CustomBusinessDay.weekmask\
        pandas.tseries.offsets.CustomBusinessHour.calendar\
        pandas.tseries.offsets.CustomBusinessHour.end\
        pandas.tseries.offsets.CustomBusinessHour.holidays\
        pandas.tseries.offsets.CustomBusinessHour.is_on_offset\
        pandas.tseries.offsets.CustomBusinessHour.n\
        pandas.tseries.offsets.CustomBusinessHour.nanos\
        pandas.tseries.offsets.CustomBusinessHour.normalize\
        pandas.tseries.offsets.CustomBusinessHour.rule_code\
        pandas.tseries.offsets.CustomBusinessHour.start\
        pandas.tseries.offsets.CustomBusinessHour.weekmask\
        pandas.tseries.offsets.CustomBusinessMonthBegin.calendar\
        pandas.tseries.offsets.CustomBusinessMonthBegin.holidays\
        pandas.tseries.offsets.CustomBusinessMonthBegin.m_offset\
        pandas.tseries.offsets.CustomBusinessMonthBegin.n\
        pandas.tseries.offsets.CustomBusinessMonthBegin.nanos\
        pandas.tseries.offsets.CustomBusinessMonthBegin.normalize\
        pandas.tseries.offsets.CustomBusinessMonthBegin.rule_code\
        pandas.tseries.offsets.CustomBusinessMonthBegin.weekmask\
        pandas.tseries.offsets.CustomBusinessMonthEnd.calendar\
        pandas.tseries.offsets.CustomBusinessMonthEnd.holidays\
        pandas.tseries.offsets.CustomBusinessMonthEnd.m_offset\
        pandas.tseries.offsets.CustomBusinessMonthEnd.n\
        pandas.tseries.offsets.CustomBusinessMonthEnd.nanos\
        pandas.tseries.offsets.CustomBusinessMonthEnd.normalize\
        pandas.tseries.offsets.CustomBusinessMonthEnd.rule_code\
        pandas.tseries.offsets.CustomBusinessMonthEnd.weekmask\
        pandas.tseries.offsets.DateOffset.is_on_offset\
        pandas.tseries.offsets.DateOffset.n\
        pandas.tseries.offsets.DateOffset.nanos\
        pandas.tseries.offsets.DateOffset.normalize\
        pandas.tseries.offsets.DateOffset.rule_code\
        pandas.tseries.offsets.Day.delta\
        pandas.tseries.offsets.Day.is_on_offset\
        pandas.tseries.offsets.Day.n\
        pandas.tseries.offsets.Day.normalize\
        pandas.tseries.offsets.Day.rule_code\
        pandas.tseries.offsets.Easter.is_on_offset\
        pandas.tseries.offsets.Easter.n\
        pandas.tseries.offsets.Easter.nanos\
        pandas.tseries.offsets.Easter.normalize\
        pandas.tseries.offsets.Easter.rule_code\
        pandas.tseries.offsets.FY5253.get_rule_code_suffix\
        pandas.tseries.offsets.FY5253.get_year_end\
        pandas.tseries.offsets.FY5253.is_anchored\
        pandas.tseries.offsets.FY5253.is_on_offset\
        pandas.tseries.offsets.FY5253.n\
        pandas.tseries.offsets.FY5253.nanos\
        pandas.tseries.offsets.FY5253.normalize\
        pandas.tseries.offsets.FY5253.rule_code\
        pandas.tseries.offsets.FY5253.startingMonth\
        pandas.tseries.offsets.FY5253.variation\
        pandas.tseries.offsets.FY5253.weekday\
        pandas.tseries.offsets.FY5253Quarter.get_rule_code_suffix\
        pandas.tseries.offsets.FY5253Quarter.get_weeks\
        pandas.tseries.offsets.FY5253Quarter.is_anchored\
        pandas.tseries.offsets.FY5253Quarter.is_on_offset\
        pandas.tseries.offsets.FY5253Quarter.n\
        pandas.tseries.offsets.FY5253Quarter.nanos\
        pandas.tseries.offsets.FY5253Quarter.normalize\
        pandas.tseries.offsets.FY5253Quarter.qtr_with_extra_week\
        pandas.tseries.offsets.FY5253Quarter.rule_code\
        pandas.tseries.offsets.FY5253Quarter.startingMonth\
        pandas.tseries.offsets.FY5253Quarter.variation\
        pandas.tseries.offsets.FY5253Quarter.weekday\
        pandas.tseries.offsets.FY5253Quarter.year_has_extra_week\
        pandas.tseries.offsets.Hour.delta\
        pandas.tseries.offsets.Hour.is_on_offset\
        pandas.tseries.offsets.Hour.n\
        pandas.tseries.offsets.Hour.normalize\
        pandas.tseries.offsets.Hour.rule_code\
        pandas.tseries.offsets.LastWeekOfMonth.is_on_offset\
        pandas.tseries.offsets.LastWeekOfMonth.n\
        pandas.tseries.offsets.LastWeekOfMonth.nanos\
        pandas.tseries.offsets.LastWeekOfMonth.normalize\
        pandas.tseries.offsets.LastWeekOfMonth.rule_code\
        pandas.tseries.offsets.LastWeekOfMonth.week\
        pandas.tseries.offsets.LastWeekOfMonth.weekday\
        pandas.tseries.offsets.Micro.delta\
        pandas.tseries.offsets.Micro.is_on_offset\
        pandas.tseries.offsets.Micro.n\
        pandas.tseries.offsets.Micro.normalize\
        pandas.tseries.offsets.Micro.rule_code\
        pandas.tseries.offsets.Milli.delta\
        pandas.tseries.offsets.Milli.is_on_offset\
        pandas.tseries.offsets.Milli.n\
        pandas.tseries.offsets.Milli.normalize\
        pandas.tseries.offsets.Milli.rule_code\
        pandas.tseries.offsets.Minute.delta\
        pandas.tseries.offsets.Minute.is_on_offset\
        pandas.tseries.offsets.Minute.n\
        pandas.tseries.offsets.Minute.normalize\
        pandas.tseries.offsets.Minute.rule_code\
        pandas.tseries.offsets.MonthBegin.is_on_offset\
        pandas.tseries.offsets.MonthBegin.n\
        pandas.tseries.offsets.MonthBegin.nanos\
        pandas.tseries.offsets.MonthBegin.normalize\
        pandas.tseries.offsets.MonthBegin.rule_code\
        pandas.tseries.offsets.MonthEnd.is_on_offset\
        pandas.tseries.offsets.MonthEnd.n\
        pandas.tseries.offsets.MonthEnd.nanos\
        pandas.tseries.offsets.MonthEnd.normalize\
        pandas.tseries.offsets.MonthEnd.rule_code\
        pandas.tseries.offsets.Nano.delta\
        pandas.tseries.offsets.Nano.is_on_offset\
        pandas.tseries.offsets.Nano.n\
        pandas.tseries.offsets.Nano.normalize\
        pandas.tseries.offsets.Nano.rule_code\
        pandas.tseries.offsets.QuarterBegin.is_anchored\
        pandas.tseries.offsets.QuarterBegin.is_on_offset\
        pandas.tseries.offsets.QuarterBegin.n\
        pandas.tseries.offsets.QuarterBegin.nanos\
        pandas.tseries.offsets.QuarterBegin.normalize\
        pandas.tseries.offsets.QuarterBegin.rule_code\
        pandas.tseries.offsets.QuarterBegin.startingMonth\
        pandas.tseries.offsets.QuarterEnd.is_anchored\
        pandas.tseries.offsets.QuarterEnd.is_on_offset\
        pandas.tseries.offsets.QuarterEnd.n\
        pandas.tseries.offsets.QuarterEnd.nanos\
        pandas.tseries.offsets.QuarterEnd.normalize\
        pandas.tseries.offsets.QuarterEnd.rule_code\
        pandas.tseries.offsets.QuarterEnd.startingMonth\
        pandas.tseries.offsets.Second.delta\
        pandas.tseries.offsets.Second.is_on_offset\
        pandas.tseries.offsets.Second.n\
        pandas.tseries.offsets.Second.normalize\
        pandas.tseries.offsets.Second.rule_code\
        pandas.tseries.offsets.SemiMonthBegin.day_of_month\
        pandas.tseries.offsets.SemiMonthBegin.is_on_offset\
        pandas.tseries.offsets.SemiMonthBegin.n\
        pandas.tseries.offsets.SemiMonthBegin.nanos\
        pandas.tseries.offsets.SemiMonthBegin.normalize\
        pandas.tseries.offsets.SemiMonthBegin.rule_code\
        pandas.tseries.offsets.SemiMonthEnd.day_of_month\
        pandas.tseries.offsets.SemiMonthEnd.is_on_offset\
        pandas.tseries.offsets.SemiMonthEnd.n\
        pandas.tseries.offsets.SemiMonthEnd.nanos\
        pandas.tseries.offsets.SemiMonthEnd.normalize\
        pandas.tseries.offsets.SemiMonthEnd.rule_code\
        pandas.tseries.offsets.Tick\
        pandas.tseries.offsets.Tick.delta\
        pandas.tseries.offsets.Tick.is_on_offset\
        pandas.tseries.offsets.Tick.n\
        pandas.tseries.offsets.Tick.normalize\
        pandas.tseries.offsets.Tick.rule_code\
        pandas.tseries.offsets.Week.is_anchored\
        pandas.tseries.offsets.Week.is_on_offset\
        pandas.tseries.offsets.Week.n\
        pandas.tseries.offsets.Week.nanos\
        pandas.tseries.offsets.Week.normalize\
        pandas.tseries.offsets.Week.rule_code\
        pandas.tseries.offsets.Week.weekday\
        pandas.tseries.offsets.WeekOfMonth.is_on_offset\
        pandas.tseries.offsets.WeekOfMonth.n\
        pandas.tseries.offsets.WeekOfMonth.nanos\
        pandas.tseries.offsets.WeekOfMonth.normalize\
        pandas.tseries.offsets.WeekOfMonth.rule_code\
        pandas.tseries.offsets.WeekOfMonth.week\
        pandas.tseries.offsets.WeekOfMonth.weekday\
        pandas.tseries.offsets.YearBegin.is_on_offset\
        pandas.tseries.offsets.YearBegin.month\
        pandas.tseries.offsets.YearBegin.n\
        pandas.tseries.offsets.YearBegin.nanos\
        pandas.tseries.offsets.YearBegin.normalize\
        pandas.tseries.offsets.YearBegin.rule_code\
        pandas.tseries.offsets.YearEnd.is_on_offset\
        pandas.tseries.offsets.YearEnd.month\
        pandas.tseries.offsets.YearEnd.n\
        pandas.tseries.offsets.YearEnd.nanos\
        pandas.tseries.offsets.YearEnd.normalize\
        pandas.tseries.offsets.YearEnd.rule_code # There should be no backslash in the final line, please keep this comment in the last ignored function
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
