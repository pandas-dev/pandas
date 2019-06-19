#!/bin/bash

set -e

export http_proxy=http://1.2.3.4 https_proxy=http://1.2.3.4;

if [ -n "$LOCALE_OVERRIDE" ]; then
    export LC_ALL="$LOCALE_OVERRIDE"
    export LANG="$LOCALE_OVERRIDE"
    PANDAS_LOCALE=`python -c 'import pandas; pandas.get_option("display.encoding")'`
    if [[ "$LOCALE_OVERRIDE" != "$PANDAS_LOCALE" ]]; then
        echo "pandas could not detect the locale. System locale: $LOCALE_OVERRIDE, pandas detected: $PANDAS_LOCALE"
        # TODO Not really aborting the tests until https://github.com/pandas-dev/pandas/issues/23923 is fixed
        # exit 1
    fi
fi

if [ "$COVERAGE" ]; then
    COVERAGE_FNAME="/tmp/test_coverage.xml"
    COVERAGE="-s --cov=pandas --cov-report=xml:$COVERAGE_FNAME"
fi

PYTEST_CMD="pytest -m \"$PATTERN\" -n auto --dist=loadfile -s --strict --durations=10 --junitxml=test-data.xml $TEST_ARGS $COVERAGE pandas"
echo $PYTEST_CMD
sh -c "$PYTEST_CMD"

if [[ "$COVERAGE" && $? == 0 ]]; then
    echo "bash <(curl -s https://codecov.io/bash) -Z -c -F $TYPE -f $COVERAGE_FNAME"
          bash <(curl -s https://codecov.io/bash) -Z -c -F $TYPE -f $COVERAGE_FNAME
fi
