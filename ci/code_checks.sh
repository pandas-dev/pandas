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

    MSG='Validate docstrings (EX02, EX04, GL01, GL02, GL03, GL04, GL05, GL06, GL07, GL09, GL10, PR03, PR04, PR05, PR06, PR08, PR09, PR10, RT01, RT02, RT04, RT05, SA02, SA03, SA04, SS01, SS02, SS03, SS04, SS05, SS06)' ; echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX02,EX04,GL01,GL02,GL03,GL04,GL05,GL06,GL07,GL09,GL10,PR03,PR04,PR05,PR06,PR08,PR09,PR10,RT01,RT02,RT04,RT05,SA02,SA03,SA04,SS01,SS02,SS03,SS04,SS05,SS06
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Partially validate docstrings (EX01)' ;  echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX01 --ignore_functions \
        pandas.Series.backfill \
        pandas.Series.pad \
        pandas.Series.hist \
        pandas.errors.AccessorRegistrationWarning \
        pandas.errors.AttributeConflictWarning \
        pandas.errors.DataError \
        pandas.errors.IncompatibilityWarning \
        pandas.errors.InvalidComparison \
        pandas.errors.IntCastingNaNError \
        pandas.errors.LossySetitemError \
        pandas.errors.MergeError \
        pandas.errors.NoBufferPresent \
        pandas.errors.NullFrequencyError \
        pandas.errors.NumbaUtilError \
        pandas.errors.OptionError \
        pandas.errors.OutOfBoundsDatetime \
        pandas.errors.OutOfBoundsTimedelta \
        pandas.errors.ParserError \
        pandas.errors.PerformanceWarning \
        pandas.errors.PyperclipException \
        pandas.errors.PyperclipWindowsException \
        pandas.errors.UnsortedIndexError \
        pandas.errors.UnsupportedFunctionCall \
        pandas.test \
        pandas.NaT \
        pandas.io.formats.style.Styler.to_html \
        pandas.read_feather \
        pandas.DataFrame.to_feather \
        pandas.read_parquet \
        pandas.read_orc \
        pandas.read_sas \
        pandas.read_spss \
        pandas.read_sql_query \
        pandas.read_gbq \
        pandas.io.stata.StataReader.data_label \
        pandas.io.stata.StataReader.value_labels \
        pandas.io.stata.StataReader.variable_labels \
        pandas.io.stata.StataWriter.write_file \
        pandas.plotting.deregister_matplotlib_converters \
        pandas.plotting.plot_params \
        pandas.plotting.register_matplotlib_converters \
        pandas.plotting.table \
        pandas.util.hash_array \
        pandas.util.hash_pandas_object \
        pandas_object \
        pandas.api.interchange.from_dataframe \
        pandas.DatetimeIndex.snap \
        pandas.core.window.rolling.Rolling.max \
        pandas.core.window.rolling.Rolling.cov \
        pandas.core.window.rolling.Rolling.skew \
        pandas.core.window.rolling.Rolling.apply \
        pandas.core.window.rolling.Window.mean \
        pandas.core.window.rolling.Window.sum \
        pandas.core.window.rolling.Window.var \
        pandas.core.window.rolling.Window.std \
        pandas.core.window.expanding.Expanding.count \
        pandas.core.window.expanding.Expanding.sum \
        pandas.core.window.expanding.Expanding.mean \
        pandas.core.window.expanding.Expanding.median \
        pandas.core.window.expanding.Expanding.min \
        pandas.core.window.expanding.Expanding.max \
        pandas.core.window.expanding.Expanding.corr \
        pandas.core.window.expanding.Expanding.cov \
        pandas.core.window.expanding.Expanding.skew \
        pandas.core.window.expanding.Expanding.apply \
        pandas.core.window.expanding.Expanding.quantile \
        pandas.core.window.ewm.ExponentialMovingWindow.mean \
        pandas.core.window.ewm.ExponentialMovingWindow.sum \
        pandas.core.window.ewm.ExponentialMovingWindow.std \
        pandas.core.window.ewm.ExponentialMovingWindow.var \
        pandas.core.window.ewm.ExponentialMovingWindow.corr \
        pandas.core.window.ewm.ExponentialMovingWindow.cov \
        pandas.api.indexers.BaseIndexer \
        pandas.api.indexers.VariableOffsetWindowIndexer \
        pandas.io.formats.style.Styler \
        pandas.io.formats.style.Styler.from_custom_template \
        pandas.io.formats.style.Styler.set_caption \
        pandas.io.formats.style.Styler.set_sticky \
        pandas.io.formats.style.Styler.set_uuid \
        pandas.io.formats.style.Styler.clear \
        pandas.io.formats.style.Styler.highlight_null \
        pandas.io.formats.style.Styler.highlight_max \
        pandas.io.formats.style.Styler.highlight_min \
        pandas.io.formats.style.Styler.bar \
        pandas.io.formats.style.Styler.to_string \
        pandas.api.extensions.ExtensionDtype \
        pandas.api.extensions.ExtensionArray \
        pandas.arrays.PandasArray \
        pandas.api.extensions.ExtensionArray._accumulate \
        pandas.api.extensions.ExtensionArray._concat_same_type \
        pandas.api.extensions.ExtensionArray._formatter \
        pandas.api.extensions.ExtensionArray._from_factorized \
        pandas.api.extensions.ExtensionArray._from_sequence \
        pandas.api.extensions.ExtensionArray._from_sequence_of_strings \
        pandas.api.extensions.ExtensionArray._hash_pandas_object \
        pandas.api.extensions.ExtensionArray._reduce \
        pandas.api.extensions.ExtensionArray._values_for_argsort \
        pandas.api.extensions.ExtensionArray._values_for_factorize \
        pandas.api.extensions.ExtensionArray.argsort \
        pandas.api.extensions.ExtensionArray.astype \
        pandas.api.extensions.ExtensionArray.copy \
        pandas.api.extensions.ExtensionArray.view \
        pandas.api.extensions.ExtensionArray.dropna \
        pandas.api.extensions.ExtensionArray.equals \
        pandas.api.extensions.ExtensionArray.factorize \
        pandas.api.extensions.ExtensionArray.fillna \
        pandas.api.extensions.ExtensionArray.insert \
        pandas.api.extensions.ExtensionArray.interpolate \
        pandas.api.extensions.ExtensionArray.isin \
        pandas.api.extensions.ExtensionArray.isna \
        pandas.api.extensions.ExtensionArray.ravel \
        pandas.api.extensions.ExtensionArray.searchsorted \
        pandas.api.extensions.ExtensionArray.shift \
        pandas.api.extensions.ExtensionArray.unique \
        pandas.api.extensions.ExtensionArray.dtype \
        pandas.api.extensions.ExtensionArray.nbytes \
        pandas.api.extensions.ExtensionArray.ndim \
        pandas.api.extensions.ExtensionArray.shape \
        pandas.api.extensions.ExtensionArray.tolist \
        pandas.DataFrame.columns \
        pandas.DataFrame.backfill \
        pandas.DataFrame.ffill \
        pandas.DataFrame.pad \
        pandas.DataFrame.swapaxes \
        pandas.DataFrame.attrs \
        pandas.DataFrame.plot \
        pandas.DataFrame.to_gbq \
        pandas.DataFrame.style \
        pandas.DataFrame.__dataframe__
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
    python doc/make.py --warnings-are-errors --single pandas.Series.value_counts
    python doc/make.py --warnings-are-errors --single pandas.Series.str.split
    python doc/make.py clean
fi

exit $RET
