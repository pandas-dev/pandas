#!/bin/bash

echo "inside $0"

source activate pandas

RET=0

if [ "$LINT" ]; then

    # pandas/_libs/src is C code, so no need to search there.
    echo "Linting *.py"
    flake8 pandas --filename=*.py --exclude pandas/_libs/src
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting *.py DONE"

    echo "Linting setup.py"
    flake8 setup.py
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting setup.py DONE"

    echo "Linting *.pyx"
    flake8 pandas --filename=*.pyx --select=E501,E302,E203,E111,E114,E221,E303,E128,E231,E126,E265,E305,E301,E127,E261,E271,E129,W291,E222,E241,E123,F403
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
    for path in '*.h' 'period_helper.c' 'datetime' 'parser' 'ujson'
    do
        echo "linting -> pandas/_libs/src/$path"
        cpplint --quiet --extensions=c,h --headers=h --filter=-readability/casting,-runtime/int,-build/include_subdir --recursive pandas/_libs/src/$path
        if [ $? -ne "0" ]; then
            RET=1
        fi
    done
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
    echo "Check for invalid testing DONE"

    echo "Check for use of lists instead of generators in built-in Python functions"

    # Example: Avoid `any([i for i in some_iterator])` in favor of `any(i for i in some_iterator)`
    #
    # Check the following functions:
    # any(), all(), sum(), max(), min(), list(), dict(), set(), frozenset(), tuple(), str.join()
    grep -R --include="*.py*" -E "[^_](any|all|sum|max|min|list|dict|set|frozenset|tuple|join)\(\[.* for .* in .*\]\)" *

    if [ $? = "0" ]; then
        RET=1
    fi
    echo "Check for use of lists instead of generators in built-in Python functions DONE"

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
else
    echo "NOT Linting"
fi

exit $RET
