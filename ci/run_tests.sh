#!/bin/bash

PATTERN=$1  # e.g. "slow and not network"
LOCALE=$2
COVERAGE=$3


# Workaround for pytest-xdist flaky collection order
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
export PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')

if [ -n "$LOCALE" ]; then
    export LC_ALL="$LOCALE"
    export LANG="$LOCALE"
    PANDAS_ENCODING=`python -c 'import pandas; pandas.get_option("display.enconding")'`
    echo "pandas detected console encoding: $PANDAS_ENCODING"
fi
if [[ "not network" == *"$PATTERN"* ]]; then
    export http_proxy=http://1.2.3.4 https_proxy=http://1.2.3.4;
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

    # TODO if $PATTERN is empty, the -m expression will probably fail
    pytest -m "$TYPE_PATTERN and $PATTERN" -n $NUM_JOBS -s --strict --durations=10 --junitxml=test-data-$TYPE.xml $COVERAGE
    if [ "$COVERAGE" && $? == 0 ]; then
        bash <(curl -s https://codecov.io/bash) -Z -c -F $TYPE -f $COVERAGE_FNAME
    fi
done
