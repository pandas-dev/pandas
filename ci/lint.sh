#!/bin/bash

echo "inside $0"

source activate pandas

RET=0

if [ "$LINT" ]; then

    # We're ignoring the following codes across the board
    #E402,  # module level import not at top of file
    #E731,  # do not assign a lambda expression, use a def
    #E741,  # do not use variables named 'l', 'O', or 'I'
    #W503,  # line break before binary operator
    #C406,  # Unnecessary (list/tuple) literal - rewrite as a dict literal.
    #C408,  # Unnecessary (dict/list/tuple) call - rewrite as a literal.
    #C409,  # Unnecessary (list/tuple) passed to tuple() - (remove the outer call to tuple()/rewrite as a tuple literal).
    #C410   # Unnecessary (list/tuple) passed to list() - (remove the outer call to list()/rewrite as a list literal).

    # pandas/_libs/src is C code, so no need to search there.
    echo "Linting *.py"
    flake8 pandas --filename=*.py --exclude pandas/_libs/src --ignore=C406,C408,C409,E402,E731,E741,W503
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting *.py DONE"

    echo "Linting setup.py"
    flake8 setup.py --ignore=E402,E731,E741,W503
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting setup.py DONE"

    echo "Linting asv_bench/benchmarks/"
    flake8 asv_bench/benchmarks/  --exclude=asv_bench/benchmarks/*.py --ignore=F811,C406,C408,C409,C410
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting asv_bench/benchmarks/*.py DONE"

    echo "Linting scripts/*.py"
    flake8 scripts --filename=*.py --ignore=C408,E402,E731,E741,W503
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting scripts/*.py DONE"

    echo "Linting doc scripts"
    flake8 doc/make.py doc/source/conf.py --ignore=E402,E731,E741,W503
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting doc scripts DONE"

    echo "Linting *.pyx"
    flake8 pandas --filename=*.pyx --select=E501,E302,E203,E111,E114,E221,E303,E128,E231,E126,E265,E305,E301,E127,E261,E271,E129,W291,E222,E241,E123,F403,C400,C401,C402,C403,C404,C407,C411
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting *.pyx DONE"

    echo "Linting *.pxi.in"
    for path in 'src'
    do
        echo "linting -> pandas/$path"
        flake8 pandas/$path --filename=*.pxi.in --select=E501,E302,E203,E111,E114,E221,E303,E231,E126,F403
        if [ $? -ne "0" ]; then
            RET=1
        fi
    done
    echo "Linting *.pxi.in DONE"

    echo "Linting *.pxd"
    for path in '_libs'
    do
        echo "linting -> pandas/$path"
        flake8 pandas/$path --filename=*.pxd --select=E501,E302,E203,E111,E114,E221,E303,E231,E126,F403
        if [ $? -ne "0" ]; then
            RET=1
        fi
    done
    echo "Linting *.pxd DONE"

    # readability/casting: Warnings about C casting instead of C++ casting
    # runtime/int: Warnings about using C number types instead of C++ ones
    # build/include_subdir: Warnings about prefacing included header files with directory

    # We don't lint all C files because we don't want to lint any that are built
    # from Cython files nor do we want to lint C files that we didn't modify for
    # this particular codebase (e.g. src/headers, src/klib, src/msgpack). However,
    # we can lint all header files since they aren't "generated" like C files are.
    echo "Linting *.c and *.h"
    for path in '*.h' 'parser' 'ujson'
    do
        echo "linting -> pandas/_libs/src/$path"
        cpplint --quiet --extensions=c,h --headers=h --filter=-readability/casting,-runtime/int,-build/include_subdir --recursive pandas/_libs/src/$path
        if [ $? -ne "0" ]; then
            RET=1
        fi
    done
    echo "linting -> pandas/_libs/tslibs/src/datetime"
    cpplint --quiet --extensions=c,h --headers=h --filter=-readability/casting,-runtime/int,-build/include_subdir --recursive pandas/_libs/tslibs/src/datetime
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting *.c and *.h DONE"

    echo "Check for invalid testing"

    # Check for the following code in testing:
    #
    # np.testing
    # np.array_equal
    grep -r -E --include '*.py' --exclude testing.py '(numpy|np)(\.testing|\.array_equal)' pandas/tests/

    if [ $? = "0" ]; then
        RET=1
    fi

    # Check for pytest.warns
    grep -r -E --include '*.py' 'pytest\.warns' pandas/tests/

    if [ $? = "0" ]; then
        RET=1
    fi

    # Check for the following code in the extension array base tests
    # tm.assert_frame_equal
    # tm.assert_series_equal
    grep -r -E --include '*.py' --exclude base.py 'tm.assert_(series|frame)_equal' pandas/tests/extension/base

    if [ $? = "0" ]; then
        RET=1
    fi

    echo "Check for invalid testing DONE"

    # Check for imports from pandas.core.common instead
    # of `import pandas.core.common as com`
    echo "Check for non-standard imports"
    grep -R --include="*.py*" -E "from pandas.core.common import " pandas
    if [ $? = "0" ]; then
        RET=1
    fi
    echo "Check for non-standard imports DONE"

    echo "Check for incorrect sphinx directives"
    SPHINX_DIRECTIVES=$(echo \
       "autosummary|contents|currentmodule|deprecated|function|image|"\
       "important|include|ipython|literalinclude|math|module|note|raw|"\
       "seealso|toctree|versionadded|versionchanged|warning" | tr -d "[:space:]")
    for path in './pandas' './doc/source'
    do
        grep -R --include="*.py" --include="*.pyx" --include="*.rst" -E "\.\. ($SPHINX_DIRECTIVES):[^:]" $path
        if [ $? = "0" ]; then
            RET=1
        fi
    done
    echo "Check for incorrect sphinx directives DONE"

    echo "Check for deprecated messages without sphinx directive"
    grep -R --include="*.py" --include="*.pyx" -E "(DEPRECATED|DEPRECATE|Deprecated)(:|,|\.)" pandas

    if [ $? = "0" ]; then
        RET=1
    fi
    echo "Check for deprecated messages without sphinx directive DONE"

    echo "Check for old-style classes"
    grep -R --include="*.py" -E "class\s\S*[^)]:" pandas scripts

    if [ $? = "0" ]; then
        RET=1
    fi
    echo "Check for old-style classes DONE"
    
    echo "Check for backticks incorrectly rendering because of missing spaces"
    grep -R --include="*.rst" -E "[a-zA-Z0-9]\`\`?[a-zA-Z0-9]" doc/source/

    if [ $? = "0" ]; then
        RET=1
    fi
    echo "Check for backticks incorrectly rendering because of missing spaces DONE"

else
    echo "NOT Linting"
fi

exit $RET
