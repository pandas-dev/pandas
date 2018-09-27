#!/bin/bash

echo "inside $0"
[[ $LINT ]] || { echo "NOT Linting"; exit 0; }

source activate pandas
RET=0


### LINTING ###

# We're ignoring the following codes across the board
# E402 module level import not at top of file
# E731 do not assign a lambda expression, use a def
# E741 do not use variables named 'l', 'O', or 'I'
# W503 line break before binary operator
# C406 Unnecessary (list/tuple) literal - rewrite as a dict literal.
# C408 Unnecessary (dict/list/tuple) call - rewrite as a literal.
# C409 Unnecessary (list/tuple) passed to tuple() - (remove the outer call to tuple()/rewrite as a tuple literal).
# C410 Unnecessary (list/tuple) passed to list() - (remove the outer call to list()/rewrite as a list literal).

# pandas/_libs/src is C code, so no need to search there.
MSG='Linting .py code' ; echo $MSG
flake8 pandas --filename=*.py --exclude pandas/_libs/src --ignore=C406,C408,C409,E402,E731,E741,W503
RET=$(($RET + $?)) ; echo $MSG "DONE"

MSG='Linting .pyx code' ; echo $MSG
flake8 pandas --filename=*.pyx --select=E501,E302,E203,E111,E114,E221,E303,E128,E231,E126,E265,E305,E301,E127,E261,E271,E129,W291,E222,E241,E123,F403,C400,C401,C402,C403,C404,C405,C406,C407,C408,C409,C410,C411
RET=$(($RET + $?)) ; echo $MSG "DONE"

MSG='Linting .pxd and .pxi.in' ; echo $MSG
flake8 pandas/_libs --filename=*.pxi.in,*.pxd --select=E501,E302,E203,E111,E114,E221,E303,E231,E126,F403
RET=$(($RET + $?)) ; echo $MSG "DONE"

MSG='Linting setup.py' ; echo $MSG
flake8 setup.py --ignore=E402,E731,E741,W503
RET=$(($RET + $?)) ; echo $MSG "DONE"

MSG='Linting scripts' ; echo $MSG
flake8 scripts --filename=*.py --ignore=C408,E402,E731,E741,W503
RET=$(($RET + $?)) ; echo $MSG "DONE"

MSG='Linting asv benchmarks' ; echo $MSG
flake8 asv_bench/benchmarks/ --exclude=asv_bench/benchmarks/*.py --ignore=F811,C406,C408,C409,C410
RET=$(($RET + $?)) ; echo $MSG "DONE"

MSG='Linting doc scripts' ; echo $MSG
flake8 doc/make.py doc/source/conf.py --ignore=E402,E731,E741,W503
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
