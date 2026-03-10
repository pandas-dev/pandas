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
    python "$BASE_DIR"/scripts/validate_docstrings.py \
        --format=actions \
        -i "pandas.IntervalDtype.subtype ES01" \
        -i "pandas.api.extensions.ExtensionArray.count ES01" \
        -i "pandas.read_html ES01" \
        -i "pandas.read_xml ES01" \
        -i "pandas.IndexSlice ES01" \
        -i "pandas.TimedeltaIndex.to_pytimedelta ES01" \
        -i "pandas.api.typing.Expanding.count ES01" \
        -i "pandas.api.typing.Expanding.sum ES01" \
        -i "pandas.api.typing.Expanding.mean ES01" \
        -i "pandas.api.typing.Expanding.median ES01" \
        -i "pandas.api.typing.Expanding.var ES01" \
        -i "pandas.api.typing.Expanding.std ES01" \
        -i "pandas.api.typing.Expanding.min ES01" \
        -i "pandas.api.typing.Expanding.max ES01" \
        -i "pandas.api.typing.Expanding.first ES01" \
        -i "pandas.api.typing.Expanding.last ES01" \
        -i "pandas.api.indexers.BaseIndexer ES01" \
        -i "pandas.api.indexers.FixedForwardWindowIndexer ES01" \
        -i "pandas.api.indexers.VariableOffsetWindowIndexer ES01" \
        -i "pandas.NamedAgg ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.idxmax ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.idxmin ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.max ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.mean ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.min ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.nunique ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.pct_change ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.prod ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.quantile ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.rank ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.rolling ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.size ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.kurt ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.sum ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.value_counts ES01" \
        -i "pandas.api.typing.SeriesGroupBy.corr ES01" \
        -i "pandas.api.typing.SeriesGroupBy.cov ES01" \
        -i "pandas.api.typing.SeriesGroupBy.is_monotonic_increasing ES01" \
        -i "pandas.api.typing.SeriesGroupBy.is_monotonic_decreasing ES01" \
        -i "pandas.api.typing.SeriesGroupBy.max ES01" \
        -i "pandas.api.typing.SeriesGroupBy.mean ES01" \
        -i "pandas.api.typing.SeriesGroupBy.min ES01" \
        -i "pandas.api.typing.SeriesGroupBy.nlargest ES01" \
        -i "pandas.api.typing.SeriesGroupBy.nsmallest ES01" \
        -i "pandas.api.typing.SeriesGroupBy.nunique ES01" \
        -i "pandas.api.typing.SeriesGroupBy.pct_change ES01" \
        -i "pandas.api.typing.SeriesGroupBy.prod ES01" \
        -i "pandas.api.typing.SeriesGroupBy.quantile ES01" \
        -i "pandas.api.typing.SeriesGroupBy.rank ES01" \
        -i "pandas.api.typing.SeriesGroupBy.rolling ES01" \
        -i "pandas.api.typing.SeriesGroupBy.size ES01" \
        -i "pandas.api.typing.SeriesGroupBy.kurt ES01" \
        -i "pandas.api.typing.SeriesGroupBy.sum ES01" \
        -i "pandas.api.typing.SeriesGroupBy.value_counts ES01" \
        -i "pandas.api.typing.DataFrameGroupBy.boxplot ES01" \
        -i "pandas.api.typing.SeriesGroupBy.hist ES01" \
        -i "pandas.Timestamp.fromisoformat SS02,SS03,ES01,SA01,EX01" \
        -i "pandas.Timedelta.resolution_string SA01" \
        -i "pandas.DatetimeIndex.asi8 GL08" \
        -i "pandas.PeriodIndex.asi8 GL08" \
        -i "pandas.TimedeltaIndex.asi8 GL08" \
        -i "pandas.DatetimeIndex.unit GL08" \
        -i "pandas.TimedeltaIndex.unit GL08" # no backslash in the last line

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
