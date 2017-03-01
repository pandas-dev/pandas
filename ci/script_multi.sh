#!/bin/bash

echo "[script multi]"

source activate pandas

# don't run the tests for the doc build
if [ x"$DOC_BUILD" != x"" ]; then
    exit 0
fi

if [ -n "$LOCALE_OVERRIDE" ]; then
    export LC_ALL="$LOCALE_OVERRIDE";
    echo "Setting LC_ALL to $LOCALE_OVERRIDE"

    pycmd='import pandas; print("pandas detected console encoding: %s" % pandas.get_option("display.encoding"))'
    python -c "$pycmd"
fi

# Workaround for pytest-xdist flaky collection order
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
export PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')
echo PYTHONHASHSEED=$PYTHONHASHSEED

if [ "$BUILD_TEST" ]; then
    echo "We are not running pytest as this is simply a build test."
elif [ "$COVERAGE" ]; then
    echo pytest -s -n 2 -m "not single" --cov=pandas --cov-append --cov-report xml:/tmp/cov.xml --junitxml=/tmp/multiple.xml $TEST_ARGS pandas
    pytest -s -n 2 -m "not single" --cov=pandas --cov-append --cov-report xml:/tmp/cov.xml --junitxml=/tmp/multiple.xml $TEST_ARGS pandas
else
    echo pytest -n 2 -m "not single" --junitxml=/tmp/multiple.xml $TEST_ARGS pandas
    pytest -n 2 -m "not single" --junitxml=/tmp/multiple.xml $TEST_ARGS pandas # TODO: doctest
fi

RET="$?"

exit "$RET"
