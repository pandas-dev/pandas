#!/bin/bash

echo "inside $0"

source activate pandas

RET=0

if [ "$LINT" ]; then
    # pandas/rpy is deprecated and will be removed.
    # pandas/src is C code, so no need to search there.
    echo "Linting  *.py"
    flake8 pandas --filename=*.py --exclude pandas/rpy,pandas/src
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting *.py DONE"

    echo "Linting *.pyx"
    flake8 pandas --filename=*.pyx --select=E501,E302,E203,E111,E114,E221,E303,E128,E231,E126
    if [ $? -ne "0" ]; then
        RET=1
    fi
    echo "Linting *.pyx DONE"

    echo "Linting *.pxi.in"
    for path in 'src'
    do
        echo "linting -> pandas/$path"
        flake8 pandas/$path --filename=*.pxi.in --select=E501,E302,E203,E111,E114,E221,E303,E231,E126
        if [ $? -ne "0" ]; then
            RET=1
        fi

    done
    echo "Linting *.pxi.in DONE"

    echo "Check for invalid testing"
    grep -r -E --include '*.py' --exclude nosetester.py --exclude testing.py '(numpy|np)\.testing' pandas
    if [ $? = "0" ]; then
        RET=1
    fi
    echo "Check for invalid testing DONE"

else
    echo "NOT Linting"
fi

exit $RET
