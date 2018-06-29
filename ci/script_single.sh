#!/bin/bash

echo "[script_single]"

source activate pandas

if [ -n "$LOCALE_OVERRIDE" ]; then
    export LC_ALL="$LOCALE_OVERRIDE";
    echo "Setting LC_ALL to $LOCALE_OVERRIDE"

    pycmd='import pandas; print("pandas detected console encoding: %s" % pandas.get_option("display.encoding"))'
    python -c "$pycmd"
fi

if [ "$SLOW" ]; then
    TEST_ARGS="--only-slow --skip-network"
fi

# Enforce absent network during testing by faking a proxy
if echo "$TEST_ARGS" | grep -e --skip-network -q; then
    export http_proxy=http://1.2.3.4 https_proxy=http://1.2.3.4;
fi

if [ "$DOC" ]; then
    echo "We are not running pytest as this is a doc-build"

elif [ "$COVERAGE" ]; then
    echo pytest -s -m "single" -r xXs --strict --cov=pandas --cov-report xml:/tmp/cov-single.xml --junitxml=/tmp/single.xml $TEST_ARGS pandas
    pytest      -s -m "single" -r xXs --strict --cov=pandas --cov-report xml:/tmp/cov-single.xml --junitxml=/tmp/single.xml $TEST_ARGS pandas

else
    echo pytest -m "single" -r xXs --junitxml=/tmp/single.xml --strict $TEST_ARGS pandas
    pytest      -m "single" -r xXs --junitxml=/tmp/single.xml --strict $TEST_ARGS pandas # TODO: doctest

fi

RET="$?"

exit "$RET"
