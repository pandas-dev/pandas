#!/bin/bash

set -e

# Workaround for pytest-xdist (it collects different tests in the workers if PYTHONHASHSEED is not set)
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
export PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')
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

time_test_file () {
    echo -n "$1 : "
    { time python -m pytest -m \"$PATTERN\" -n 1 -s --strict --durations=10 --junitxml=test-data.xml $TEST_ARGS $1 > /dev/null ; } 2>&1 | grep "real" | cut -f2
}
export -f time_test_file
find pandas -name "test_*.py" -exec bash -c 'time_test_file "$0"' {} \;
