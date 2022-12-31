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
#   $ ./ci/code_checks.sh lint          # run linting only
#   $ ./ci/code_checks.sh patterns      # check for patterns that should not exist
#   $ ./ci/code_checks.sh code          # checks on imported code
#   $ ./ci/code_checks.sh doctests      # run doctests
#   $ ./ci/code_checks.sh docstrings    # validate docstring errors
#   $ ./ci/code_checks.sh single-docs   # check single-page docs build warning-free
#   $ ./ci/code_checks.sh notebooks     # check execution of documentation notebooks

[[ -z "$1" || "$1" == "lint" || "$1" == "patterns" || "$1" == "code" || "$1" == "doctests" || "$1" == "docstrings" || "$1" == "single-docs" || "$1" == "notebooks" ]] || \
    { echo "Unknown command $1. Usage: $0 [lint|patterns|code|doctests|docstrings|single-docs|notebooks]"; exit 9999; }

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

    MSG='Check for use of exec' ; echo $MSG
    invgrep -R --include="*.py*" -E "[^a-zA-Z0-9_]exec\(" pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for pytest warns' ; echo $MSG
    invgrep -r -E --include '*.py' 'pytest\.warns' pandas/tests/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for pytest raises without context' ; echo $MSG
    invgrep -r -E --include '*.py' "[[:space:]] pytest.raises" pandas/tests/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for use of builtin filter function' ; echo $MSG
    invgrep -R --include="*.py" -P '(?<!def)[\(\s]filter\(' pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    # Check for the following code in testing: `np.testing` and `np.array_equal`
    MSG='Check for invalid testing' ; echo $MSG
    invgrep -r -E --include '*.py' --exclude testing.py '(numpy|np)(\.testing|\.array_equal)' pandas/tests/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

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

    # Check for the following code in testing: `unittest.mock`, `mock.Mock()` or `mock.patch`
    MSG='Check that unittest.mock is not used (pytest builtin monkeypatch fixture should be used instead)' ; echo $MSG
    invgrep -r -E --include '*.py' '(unittest(\.| import )mock|mock\.Mock\(\)|mock\.patch)' pandas/tests/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for use of {foo!r} instead of {repr(foo)}' ; echo $MSG
    invgrep -R --include=*.{py,pyx} '!r}' pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"
    echo $MSG "DONE"
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

    MSG='Validate docstrings (EX04, GL01, GL02, GL03, GL04, GL05, GL06, GL07, GL09, GL10, PR03, PR04, PR05, PR06, PR08, PR09, PR10, RT01, RT04, RT05, SA02, SA03, SA04, SS01, SS02, SS03, SS04, SS05, SS06)' ; echo $MSG
    $BASE_DIR/scripts/validate_docstrings.py --format=actions --errors=EX04,GL01,GL02,GL03,GL04,GL05,GL06,GL07,GL09,GL10,PR03,PR04,PR05,PR06,PR08,PR09,PR10,RT01,RT04,RT05,SA02,SA03,SA04,SS01,SS02,SS03,SS04,SS05,SS06
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
