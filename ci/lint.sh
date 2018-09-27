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
