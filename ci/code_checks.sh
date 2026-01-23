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
    { echo "Unknown command $1. Usage: $0 [code|doctests|docstrings|single-docs|notebooks]"; exit 1; }

BASE_DIR="$(dirname "$0")/.."
RET=0

### CODE ###
if [[ -z "$CHECK" || "$CHECK" == "code" ]]; then

    MSG='Check import. No warnings, and blocklist some optional dependencies' ; echo "$MSG"
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
    RET=$(($RET + $?)) ; echo "$MSG" "DONE"

fi

### DOCTESTS ###
if [[ -z "$CHECK" || "$CHECK" == "doctests" ]]; then

    MSG='Python and Cython Doctests' ; echo "$MSG"
    # Using future.python_scalars=True avoids NumPy scalar reprs in docstrings.
    PANDAS_FUTURE_PYTHON_SCALARS="1" python -c 'import pandas as pd; pd.test(run_doctests=True)'
    RET=$(($RET + $?)) ; echo "$MSG" "DONE"

fi

### DOCSTRINGS ###
if [[ -z "$CHECK" || "$CHECK" == "docstrings" ]]; then

    MSG='Validate Docstrings' ; echo "$MSG"
    "$BASE_DIR"/scripts/validate_docstrings.py \
        --format=actions \
        -i ES01 `# For now it is ok if docstrings are missing the extended summary` \
        -i "pandas.Series.dt PR01" `# Accessors are implemented as classes, but we do not document the Parameters section` \
        -i "pandas.tseries.offsets.BDay PR02,SA01" \
        -i "pandas.tseries.offsets.BusinessDay PR02,SA01" \
        -i "pandas.tseries.offsets.BusinessHour PR02,SA01" \
        -i "pandas.tseries.offsets.CBMonthBegin PR02" \
        -i "pandas.tseries.offsets.CBMonthEnd PR02" \
        -i "pandas.tseries.offsets.CDay PR02,SA01" \
        -i "pandas.tseries.offsets.CustomBusinessDay PR02,SA01" \
        -i "pandas.tseries.offsets.CustomBusinessHour PR02,SA01" \
        -i "pandas.tseries.offsets.CustomBusinessMonthBegin PR02" \
        -i "pandas.tseries.offsets.CustomBusinessMonthEnd PR02" \
        -i "pandas.tseries.offsets.Day.is_on_offset GL08" \
        -i "pandas.tseries.offsets.Easter.is_on_offset GL08" \
        -i "pandas.tseries.offsets.FY5253.is_on_offset GL08" \
        -i "pandas.tseries.offsets.FY5253.startingMonth GL08" \
        -i "pandas.tseries.offsets.FY5253Quarter.get_weeks GL08" \
        -i "pandas.tseries.offsets.FY5253Quarter.is_on_offset GL08" \
        -i "pandas.tseries.offsets.FY5253Quarter.rule_code GL08" \
        -i "pandas.tseries.offsets.FY5253Quarter.startingMonth GL08" \
        -i "pandas.tseries.offsets.LastWeekOfMonth.is_on_offset GL08" \
        -i "pandas.tseries.offsets.LastWeekOfMonth.week GL08" \
        -i "pandas.tseries.offsets.LastWeekOfMonth.weekday GL08" \
        -i "pandas.tseries.offsets.SemiMonthBegin.is_on_offset GL08" \
        -i "pandas.tseries.offsets.SemiMonthBegin.rule_code GL08" \
        -i "pandas.tseries.offsets.SemiMonthEnd.is_on_offset GL08" \
        -i "pandas.tseries.offsets.SemiMonthEnd.rule_code GL08" \
        -i "pandas.tseries.offsets.Week.is_on_offset GL08" \
        -i "pandas.tseries.offsets.WeekOfMonth.is_on_offset GL08" \
        -i "pandas.tseries.offsets.WeekOfMonth.week GL08" \
        -i "pandas.tseries.offsets.WeekOfMonth.weekday GL08" # There should be no backslash in the final line, please keep this comment in the last ignored function

    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCUMENTATION NOTEBOOKS ###
if [[ -z "$CHECK" || "$CHECK" == "notebooks" ]]; then

    MSG='Notebooks' ; echo $MSG
    jupyter nbconvert --execute "$(find doc/source -name '*.ipynb')" --to notebook
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### SINGLE-PAGE DOCS ###
if [[ -z "$CHECK" || "$CHECK" == "single-docs" ]]; then
    python doc/make.py --warnings-are-errors --no-browser --single pandas.Series.value_counts
    python doc/make.py --warnings-are-errors --no-browser --single pandas.Series.str.split
fi

exit $RET
