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

    MSG='Doctests' ; echo $MSG
    # Ignore test_*.py files or else the unit tests will run
    python -m pytest --doctest-modules --ignore-glob="**/test_*.py" pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Cython Doctests' ; echo $MSG
    python -m pytest --doctest-cython pandas/_libs
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCSTRINGS ###
if [[ -z "$CHECK" || "$CHECK" == "docstrings" ]]; then

    MSG='Validate docstrings (EX04, GL01, GL02, GL03, GL04, GL05, GL06, GL07, GL09, GL10, PR03, PR04, PR05, PR06, PR08, PR09, PR10, RT01, RT02, RT04, RT05, SA02, SA03, SA04, SS01, SS02, SS03, SS04, SS05, SS06)' ; echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX04,GL01,GL02,GL03,GL04,GL05,GL06,GL07,GL09,GL10,PR03,PR04,PR05,PR06,PR08,PR09,PR10,RT01,RT02,RT04,RT05,SA02,SA03,SA04,SS01,SS02,SS03,SS04,SS05,SS06
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
