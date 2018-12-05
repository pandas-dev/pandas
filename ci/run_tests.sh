#!/bin/bash

set -e

if [ "$DOC" ]; then
    echo "We are not running pytest as this is a doc-build"
    exit 0
fi

# Workaround for pytest-xdist flaky collection order
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
export PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')

if [ -n "$LOCALE_OVERRIDE" ]; then
    export LC_ALL="$LOCALE_OVERRIDE"
    export LANG="$LOCALE_OVERRIDE"
    PANDAS_LOCALE=`python -c 'import pandas; pandas.get_option("display.encoding")'`
    if [[ "$LOCALE_OVERIDE" != "$PANDAS_LOCALE" ]]; then
        echo "pandas could not detect the locale. System locale: $LOCALE_OVERRIDE, pandas detected: $PANDAS_LOCALE"
        # TODO Not really aborting the tests until https://github.com/pandas-dev/pandas/issues/23923 is fixed
        # exit 1
    fi
fi
if [[ "not network" == *"$PATTERN"* ]]; then
    export http_proxy=http://1.2.3.4 https_proxy=http://1.2.3.4;
fi

PYTEST="pytest -m \"$PATTERN\" --junitxml=test-data.xml $TEST_ARGS"
if [ "$COVERAGE" ]; then
    COVERAGE_FNAME="/tmp/coverage.xml"
    $PYTEST --cov=pandas --cov-report=xml:$COVERAGE_FNAME pandas
    bash <(curl -s https://codecov.io/bash) -Z -c -f $COVERAGE_FNAME
else
    $PYTEST pandas
fi
