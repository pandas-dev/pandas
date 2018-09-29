#!/bin/bash

echo "inside $0"
[[ $LINT ]] || { echo "NOT Linting"; exit 0; }

source activate pandas
RET=0


### LINTING ###

# `setup.cfg` contains the list of error codes that are being ignored in flake8

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


### CHECKS ###

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


exit $RET
