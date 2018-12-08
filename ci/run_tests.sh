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


if [ -n "$PATTERN" ]; then
    PATTERN=" and $PATTERN"
fi

for TYPE in single multiple
do
    if [ "$COVERAGE" ]; then
        COVERAGE_FNAME="/tmp/coc-$TYPE.xml"
        COVERAGE="-s --cov=pandas --cov-report=xml:$COVERAGE_FNAME"
    fi

    TYPE_PATTERN=$TYPE
    NUM_JOBS=1
    if [[ "$TYPE_PATTERN" == "multiple" ]]; then
        TYPE_PATTERN="not single"
        NUM_JOBS=2
    fi

    PYTEST_CMD="pytest -m \"$TYPE_PATTERN$PATTERN\" -n $NUM_JOBS -s --strict --durations=10 --junitxml=test-data-$TYPE.xml $TEST_ARGS $COVERAGE pandas"
    echo $PYTEST_CMD
    # if no tests are found (the case of "single and slow"), pytest exits with code 5, and would make the script fail, if not for the below code
    sh -c "$PYTEST_CMD; ret=\$?; [ \$ret = 5 ] && exit 0 || exit \$ret"

    if [[ "$COVERAGE" && $? == 0 ]]; then
        echo "uploading coverage for $TYPE tests"
        echo "bash <(curl -s https://codecov.io/bash) -Z -c -F $TYPE -f $COVERAGE_FNAME"
              bash <(curl -s https://codecov.io/bash) -Z -c -F $TYPE -f $COVERAGE_FNAME
    fi
done
