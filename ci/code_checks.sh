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
        pandas.Series.plot.bar\
        pandas.Series.plot.barh\
        pandas.Series.plot.line\
        pandas.Series.plot.pie\
        pandas.DataFrame.plot\
        pandas.DataFrame.plot.bar\
        pandas.DataFrame.plot.barh\
        pandas.DataFrame.plot.line\
        pandas.DataFrame.plot.pie\
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
        pandas.core.groupby.SeriesGroupBy.apply\
        pandas.core.groupby.DataFrameGroupBy.nth\
        pandas.core.groupby.DataFrameGroupBy.rolling\
        pandas.core.groupby.SeriesGroupBy.nth\
        pandas.core.groupby.SeriesGroupBy.rolling\
        pandas.core.groupby.DataFrameGroupBy.hist\
        pandas.core.groupby.DataFrameGroupBy.plot\
        pandas.core.groupby.SeriesGroupBy.plot # There should be no backslash in the final line, please keep this comment in the last ignored function
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
