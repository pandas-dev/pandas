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

[[ -z "$1" || "$1" == "code" || "$1" == "doctests" || "$1" == "docstrings" || "$1" == "single-docs" || "$1" == "notebooks" ]] || \
    { echo "Unknown command $1. Usage: $0 [code|doctests|docstrings|single-docs|notebooks]"; exit 9999; }

BASE_DIR="$(dirname $0)/.."
RET=0
CHECK=$1

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

    MSG='Validate docstrings (EX01, EX02, EX04, GL01, GL02, GL03, GL04, GL05, GL06, GL07, GL09, GL10, PR03, PR04, PR05, PR06, PR08, PR09, PR10, RT01, RT02, RT04, RT05, SA02, SA03, SA04, SS01, SS02, SS03, SS04, SS05, SS06)' ; echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX01,EX02,EX04,GL01,GL02,GL03,GL04,GL05,GL06,GL07,GL09,GL10,PR03,PR04,PR05,PR06,PR08,PR09,PR10,RT01,RT02,RT04,RT05,SA02,SA03,SA04,SS01,SS02,SS03,SS04,SS05,SS06
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Partially validate docstrings (EX03)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX03 --ignore_functions \
        pandas.Series.dt.day_name \
        pandas.Series.str.len \
        pandas.Series.cat.set_categories \
        pandas.Series.plot.bar \
        pandas.Series.plot.hist \
        pandas.Series.plot.line \
        pandas.Series.to_sql \
        pandas.Series.to_latex \
        pandas.errors.CategoricalConversionWarning \
        pandas.errors.ChainedAssignmentError \
        pandas.errors.ClosedFileError \
        pandas.errors.DatabaseError \
        pandas.errors.IndexingError \
        pandas.errors.InvalidColumnName \
        pandas.errors.NumExprClobberingError \
        pandas.errors.PossibleDataLossError \
        pandas.errors.PossiblePrecisionLoss \
        pandas.errors.SettingWithCopyError \
        pandas.errors.SettingWithCopyWarning \
        pandas.errors.SpecificationError \
        pandas.errors.UndefinedVariableError \
        pandas.errors.ValueLabelTypeMismatch \
        pandas.Timestamp.ceil \
        pandas.Timestamp.floor \
        pandas.Timestamp.round \
        pandas.read_pickle \
        pandas.ExcelWriter \
        pandas.read_json \
        pandas.io.json.build_table_schema \
        pandas.DataFrame.to_latex \
        pandas.io.formats.style.Styler.to_latex \
        pandas.read_parquet \
        pandas.DataFrame.to_sql \
        pandas.read_stata \
        pandas.core.resample.Resampler.pipe \
        pandas.core.resample.Resampler.fillna \
        pandas.core.resample.Resampler.interpolate \
        pandas.plotting.scatter_matrix \
        pandas.pivot \
        pandas.merge_asof \
        pandas.wide_to_long \
        pandas.Index.rename \
        pandas.Index.droplevel \
        pandas.Index.isin \
        pandas.CategoricalIndex.set_categories \
        pandas.MultiIndex.names \
        pandas.MultiIndex.droplevel \
        pandas.IndexSlice \
        pandas.DatetimeIndex.month_name \
        pandas.DatetimeIndex.day_name \
        pandas.core.window.rolling.Rolling.corr \
        pandas.Grouper \
        pandas.core.groupby.SeriesGroupBy.apply \
        pandas.core.groupby.DataFrameGroupBy.apply \
        pandas.core.groupby.SeriesGroupBy.transform \
        pandas.core.groupby.SeriesGroupBy.pipe \
        pandas.core.groupby.DataFrameGroupBy.pipe \
        pandas.core.groupby.DataFrameGroupBy.describe \
        pandas.core.groupby.DataFrameGroupBy.idxmax \
        pandas.core.groupby.DataFrameGroupBy.idxmin \
        pandas.core.groupby.DataFrameGroupBy.value_counts \
        pandas.core.groupby.SeriesGroupBy.describe \
        pandas.core.groupby.DataFrameGroupBy.boxplot \
        pandas.core.groupby.DataFrameGroupBy.hist \
        pandas.io.formats.style.Styler.map \
        pandas.io.formats.style.Styler.apply_index \
        pandas.io.formats.style.Styler.map_index \
        pandas.io.formats.style.Styler.format \
        pandas.io.formats.style.Styler.format_index \
        pandas.io.formats.style.Styler.relabel_index \
        pandas.io.formats.style.Styler.hide \
        pandas.io.formats.style.Styler.set_td_classes \
        pandas.io.formats.style.Styler.set_tooltips \
        pandas.io.formats.style.Styler.set_uuid \
        pandas.io.formats.style.Styler.pipe \
        pandas.io.formats.style.Styler.highlight_between \
        pandas.io.formats.style.Styler.highlight_quantile \
        pandas.io.formats.style.Styler.background_gradient \
        pandas.io.formats.style.Styler.text_gradient \
        pandas.DataFrame.values \
        pandas.DataFrame.loc \
        pandas.DataFrame.iloc \
        pandas.DataFrame.groupby \
        pandas.DataFrame.describe \
        pandas.DataFrame.skew \
        pandas.DataFrame.var \
        pandas.DataFrame.idxmax \
        pandas.DataFrame.idxmin \
        pandas.DataFrame.last \
        pandas.DataFrame.pivot \
        pandas.DataFrame.sort_values \
        pandas.DataFrame.tz_convert \
        pandas.DataFrame.tz_localize \
        pandas.DataFrame.plot.bar \
        pandas.DataFrame.plot.hexbin \
        pandas.DataFrame.plot.hist \
        pandas.DataFrame.plot.line \
        pandas.DataFrame.hist \
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
