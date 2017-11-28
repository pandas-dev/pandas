#!/bin/bash -e

echo "[script multi]"

source activate pandas

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
    echo "[build-test]"

    echo "[env]"
    pip list --format columns |grep pandas

    echo "[running]"
    cd /tmp
    unset PYTHONPATH

    echo "[build-test: single]"
    python -c 'import pandas; pandas.test(["--skip-slow", "--skip-network", "-r xX", "-m single"])'

    echo "[build-test: not single]"
    python -c 'import pandas; pandas.test(["-n 2", "--skip-slow", "--skip-network", "-r xX", "-m not single"])'

elif [ "$DOC" ]; then
    echo "We are not running pytest as this is a doc-build"

elif [ "$COVERAGE" ]; then
    echo pytest -s -n 2 -m "not single" --cov=pandas --cov-report xml:/tmp/cov-multiple.xml --junitxml=/tmp/multiple.xml --strict $TEST_ARGS pandas
    pytest -s -n 2 -m "not single" --cov=pandas --cov-report xml:/tmp/cov-multiple.xml --junitxml=/tmp/multiple.xml --strict $TEST_ARGS pandas

elif [ "$SLOW" ]; then
    TEST_ARGS="--only-slow --skip-network"
    echo pytest -r xX -m "not single and slow" -v --junitxml=/tmp/multiple.xml --strict $TEST_ARGS pandas
    pytest -r xX -m "not single and slow" -v --junitxml=/tmp/multiple.xml --strict $TEST_ARGS pandas

else
    echo pytest -n 2 -r xX -m "not single" --junitxml=/tmp/multiple.xml --strict $TEST_ARGS pandas
    pytest -n 2 -r xX -m "not single" --junitxml=/tmp/multiple.xml --strict $TEST_ARGS pandas # TODO: doctest

fi

RET="$?"

exit "$RET"
