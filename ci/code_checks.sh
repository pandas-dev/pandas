#!/bin/bash
#
# Run checks related to code quality.
#
# This script is intended for both the CI and to check locally that code standards are
# respected. We are currently linting (PEP-8 and similar), looking for patterns of
# common mistakes (sphinx directives with missing blank lines, old style classes,
# unwanted imports...), we run doctests here (currently some files only), and we
# validate formatting error in docstrings.
#
# Usage:
#   $ ./ci/code_checks.sh               # run all checks
#   $ ./ci/code_checks.sh lint          # run linting only
#   $ ./ci/code_checks.sh patterns      # check for patterns that should not exist
#   $ ./ci/code_checks.sh code          # checks on imported code
#   $ ./ci/code_checks.sh doctests      # run doctests
#   $ ./ci/code_checks.sh docstrings    # validate docstring errors
#   $ ./ci/code_checks.sh typing        # run static type analysis

[[ -z "$1" || "$1" == "lint" || "$1" == "patterns" || "$1" == "code" || "$1" == "doctests" || "$1" == "docstrings" || "$1" == "typing" ]] || \
    { echo "Unknown command $1. Usage: $0 [lint|patterns|code|doctests|docstrings|typing]"; exit 9999; }

BASE_DIR="$(dirname $0)/.."
RET=0
CHECK=$1

function invgrep {
    # grep with inverse exist status and formatting for azure-pipelines
    #
    # This function works exactly as grep, but with opposite exit status:
    # - 0 (success) when no patterns are found
    # - 1 (fail) when the patterns are found
    #
    # This is useful for the CI, as we want to fail if one of the patterns
    # that we want to avoid is found by grep.
    grep -n "$@" | sed "s/^/$INVGREP_PREPEND/" | sed "s/$/$INVGREP_APPEND/" ; EXIT_STATUS=${PIPESTATUS[0]}
    return $((! $EXIT_STATUS))
}

if [[ "$GITHUB_ACTIONS" == "true" ]]; then
    INVGREP_PREPEND="##[error]"
fi

### LINTING ###
if [[ -z "$CHECK" || "$CHECK" == "lint" ]]; then

    # Check that cython casting is of the form `<type>obj` as opposed to `<type> obj`;
    # it doesn't make a difference, but we want to be internally consistent.
    # Note: this grep pattern is (intended to be) equivalent to the python
    # regex r'(?<![ ->])> '
    MSG='Linting .pyx code for spacing conventions in casting' ; echo $MSG
    invgrep -r -E --include '*.pyx' --include '*.pxi.in' '[a-zA-Z0-9*]> ' pandas/_libs
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    # readability/casting: Warnings about C casting instead of C++ casting
    # runtime/int: Warnings about using C number types instead of C++ ones
    # build/include_subdir: Warnings about prefacing included header files with directory

fi

### PATTERNS ###
if [[ -z "$CHECK" || "$CHECK" == "patterns" ]]; then

    # Check for the following code in the extension array base tests: `tm.assert_frame_equal` and `tm.assert_series_equal`
    MSG='Check for invalid EA testing' ; echo $MSG
    invgrep -r -E --include '*.py' --exclude base.py 'tm.assert_(series|frame)_equal' pandas/tests/extension/base
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for deprecated messages without sphinx directive' ; echo $MSG
    invgrep -R --include="*.py" --include="*.pyx" -E "(DEPRECATED|DEPRECATE|Deprecated)(:|,|\.)" pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for backticks incorrectly rendering because of missing spaces' ; echo $MSG
    invgrep -R --include="*.rst" -E "[a-zA-Z0-9]\`\`?[a-zA-Z0-9]" doc/source/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for unnecessary random seeds in asv benchmarks' ; echo $MSG
    invgrep -R --exclude pandas_vb_common.py -E 'np.random.seed' asv_bench/benchmarks/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### CODE ###
if [[ -z "$CHECK" || "$CHECK" == "code" ]]; then

    MSG='Check import. No warnings, and blocklist some optional dependencies' ; echo $MSG
    python -W error -c "
import sys
import pandas

blocklist = {'bs4', 'gcsfs', 'html5lib', 'http', 'ipython', 'jinja2', 'hypothesis',
             'lxml', 'matplotlib', 'openpyxl', 'py', 'pytest', 's3fs', 'scipy',
             'tables', 'urllib.request', 'xlrd', 'xlsxwriter', 'xlwt'}

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

    MSG='Doctests for individual files' ; echo $MSG
    pytest -q --doctest-modules \
      pandas/core/accessor.py \
      pandas/core/aggregation.py \
      pandas/core/algorithms.py \
      pandas/core/base.py \
      pandas/core/construction.py \
      pandas/core/frame.py \
      pandas/core/generic.py \
      pandas/core/indexers.py \
      pandas/core/nanops.py \
      pandas/core/series.py \
      pandas/io/sql.py
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Doctests for directories' ; echo $MSG
    pytest -q --doctest-modules \
      pandas/_libs/ \
      pandas/api/ \
      pandas/arrays/ \
      pandas/compat/ \
      pandas/core/array_algos/ \
      pandas/core/arrays/ \
      pandas/core/computation/ \
      pandas/core/dtypes/ \
      pandas/core/groupby/ \
      pandas/core/indexes/ \
      pandas/core/ops/ \
      pandas/core/reshape/ \
      pandas/core/strings/ \
      pandas/core/tools/ \
      pandas/core/window/ \
      pandas/errors/ \
      pandas/io/clipboard/ \
      pandas/io/json/ \
      pandas/io/excel/ \
      pandas/io/parsers/ \
      pandas/io/sas/ \
      pandas/tseries/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCSTRINGS ###
if [[ -z "$CHECK" || "$CHECK" == "docstrings" ]]; then

    MSG='Validate docstrings (GL03, GL04, GL05, GL06, GL07, GL09, GL10, SS01, SS02, SS04, SS05, PR03, PR04, PR05, PR10, EX04, RT01, RT04, RT05, SA02, SA03)' ; echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=GL03,GL04,GL05,GL06,GL07,GL09,GL10,SS02,SS04,SS05,PR03,PR04,PR05,PR10,EX04,RT01,RT04,RT05,SA02,SA03
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### TYPING ###
if [[ -z "$CHECK" || "$CHECK" == "typing" ]]; then

    echo "mypy --version"
    mypy --version

    MSG='Performing static analysis using mypy' ; echo $MSG
    mypy pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"
fi

exit $RET
