#!/bin/bash
#
# Run checks related to code quality.
#
# This script is intended for both the CI and to check locally that code standards are
# respected. We are currently linting (PEP-8 and similar), looking for patterns of
# common mistakes (sphinx directives with missing blank lines, old style classes,
# unwanted imports...), and we also run doctests here (currently some files only).
# In the future we may want to add the validation of docstrings and other checks here.
#
# Usage:
#   $ ./ci/code_checks.sh             # run all checks
#   $ ./ci/code_checks.sh lint        # run linting only
#   $ ./ci/code_checks.sh patterns    # check for patterns that should not exist
#   $ ./ci/code_checks.sh doctests    # run doctests

echo "inside $0"
[[ $LINT ]] || { echo "NOT Linting. To lint use: LINT=true $0 $1"; exit 0; }
[[ -z "$1" || "$1" == "lint" || "$1" == "patterns" || "$1" == "doctests" ]] || { echo "Unkown command $1. Usage: $0 [lint|patterns|doctests]"; exit 9999; }

source activate pandas
RET=0
CHECK=$1


### LINTING ###
if [[ -z "$CHECK" || "$CHECK" == "lint" ]]; then

    # `setup.cfg` contains the list of error codes that are being ignored in flake8

    echo "flake8 --version"
    flake8 --version

    # pandas/_libs/src is C code, so no need to search there.
    MSG='Linting .py code' ; echo $MSG
    flake8 .
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Linting .pyx code' ; echo $MSG
    flake8 pandas --filename=*.pyx --select=E501,E302,E203,E111,E114,E221,E303,E128,E231,E126,E265,E305,E301,E127,E261,E271,E129,W291,E222,E241,E123,F403,C400,C401,C402,C403,C404,C405,C406,C407,C408,C409,C410,C411
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Linting .pxd and .pxi.in' ; echo $MSG
    flake8 pandas/_libs --filename=*.pxi.in,*.pxd --select=E501,E302,E203,E111,E114,E221,E303,E231,E126,F403
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    # readability/casting: Warnings about C casting instead of C++ casting
    # runtime/int: Warnings about using C number types instead of C++ ones
    # build/include_subdir: Warnings about prefacing included header files with directory

    # We don't lint all C files because we don't want to lint any that are built
    # from Cython files nor do we want to lint C files that we didn't modify for
    # this particular codebase (e.g. src/headers, src/klib, src/msgpack). However,
    # we can lint all header files since they aren't "generated" like C files are.
    MSG='Linting .c and .h' ; echo $MSG
    cpplint --quiet --extensions=c,h --headers=h --recursive --filter=-readability/casting,-runtime/int,-build/include_subdir pandas/_libs/src/*.h pandas/_libs/src/parser pandas/_libs/ujson pandas/_libs/tslibs/src/datetime
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    # Imports - Check formatting using isort see setup.cfg for settings
    MSG='Check import format using isort ' ; echo $MSG
    isort --recursive --check-only pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### PATTERNS ###
if [[ -z "$CHECK" || "$CHECK" == "patterns" ]]; then

    # Check for imports from pandas.core.common instead of `import pandas.core.common as com`
    MSG='Check for non-standard imports' ; echo $MSG
    ! grep -R --include="*.py*" -E "from pandas.core.common import " pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for pytest warns' ; echo $MSG
    ! grep -r -E --include '*.py' 'pytest\.warns' pandas/tests/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    # Check for the following code in testing: `np.testing` and `np.array_equal`
    MSG='Check for invalid testing' ; echo $MSG
    ! grep -r -E --include '*.py' --exclude testing.py '(numpy|np)(\.testing|\.array_equal)' pandas/tests/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    # Check for the following code in the extension array base tests: `tm.assert_frame_equal` and `tm.assert_series_equal`
    MSG='Check for invalid EA testing' ; echo $MSG
    ! grep -r -E --include '*.py' --exclude base.py 'tm.assert_(series|frame)_equal' pandas/tests/extension/base
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for deprecated messages without sphinx directive' ; echo $MSG
    ! grep -R --include="*.py" --include="*.pyx" -E "(DEPRECATED|DEPRECATE|Deprecated)(:|,|\.)" pandas
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for old-style classes' ; echo $MSG
    ! grep -R --include="*.py" -E "class\s\S*[^)]:" pandas scripts
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for backticks incorrectly rendering because of missing spaces' ; echo $MSG
    ! grep -R --include="*.rst" -E "[a-zA-Z0-9]\`\`?[a-zA-Z0-9]" doc/source/
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for incorrect sphinx directives' ; echo $MSG
    ! grep -R --include="*.py" --include="*.pyx" --include="*.rst" -E "\.\. (autosummary|contents|currentmodule|deprecated|function|image|important|include|ipython|literalinclude|math|module|note|raw|seealso|toctree|versionadded|versionchanged|warning):[^:]" ./pandas ./doc/source
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Check for modules that pandas should not import' ; echo $MSG
    python -c "
import sys
import pandas

blacklist = {'bs4', 'gcsfs', 'html5lib', 'ipython', 'jinja2' 'hypothesis',
             'lxml', 'numexpr', 'openpyxl', 'py', 'pytest', 's3fs', 'scipy',
             'tables', 'xlrd', 'xlsxwriter', 'xlwt'}
mods = blacklist & set(m.split('.')[0] for m in sys.modules)
if mods:
    sys.stderr.write('pandas should not import: {}\n'.format(', '.join(mods)))
    sys.exit(len(mods))
    "
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

### DOCTESTS ###
if [[ -z "$CHECK" || "$CHECK" == "doctests" ]]; then

    MSG='Doctests frame.py' ; echo $MSG
    pytest -q --doctest-modules pandas/core/frame.py \
        -k"-axes -combine -itertuples -join -nlargest -nsmallest -nunique -pivot_table -quantile -query -reindex -reindex_axis -replace -round -set_index -stack -to_stata"
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Doctests series.py' ; echo $MSG
    pytest -q --doctest-modules pandas/core/series.py \
        -k"-nonzero -reindex -searchsorted -to_dict"
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Doctests generic.py' ; echo $MSG
    pytest -q --doctest-modules pandas/core/generic.py \
        -k"-_set_axis_name -_xs -describe -droplevel -groupby -interpolate -pct_change -pipe -reindex -reindex_axis -resample -to_json -transpose -values -xs"
    RET=$(($RET + $?)) ; echo $MSG "DONE"

    MSG='Doctests top-level reshaping functions' ; echo $MSG
    pytest -q --doctest-modules \
        pandas/core/reshape/concat.py \
        pandas/core/reshape/pivot.py \
        pandas/core/reshape/reshape.py \
        pandas/core/reshape/tile.py \
        -k"-crosstab -pivot_table -cut"
    RET=$(($RET + $?)) ; echo $MSG "DONE"

fi

exit $RET
