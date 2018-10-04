#!/bin/bash

echo "[script_single]"

source activate pandas

if [ -n "$LOCALE_OVERRIDE" ]; then
    echo "Setting LC_ALL and LANG to $LOCALE_OVERRIDE"
    export LC_ALL="$LOCALE_OVERRIDE";
    export LANG="$LOCALE_OVERRIDE";

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
    echo pytest -s -m "single" -r xXs --strict --cov=pandas --cov-report xml:/tmp/cov-single.xml --junitxml=test-data-single.xml $TEST_ARGS pandas
    pytest      -s -m "single" -r xXs --strict --cov=pandas --cov-report xml:/tmp/cov-single.xml --junitxml=test-data-single.xml $TEST_ARGS pandas

    echo pytest -s -r xXs --strict scripts
    pytest      -s -r xXs --strict scripts
else
    echo pytest -m "single" -r xXs --junitxml=test-data-single.xml --strict $TEST_ARGS pandas
    pytest      -m "single" -r xXs --junitxml=test-data-single.xml --strict $TEST_ARGS pandas

fi

RET="$?"

exit "$RET"
