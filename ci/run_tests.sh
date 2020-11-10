#!/bin/bash -e

# Workaround for pytest-xdist (it collects different tests in the workers if PYTHONHASHSEED is not set)
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
export PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')

if [[ "not network" == *"$PATTERN"* ]]; then
    export http_proxy=http://1.2.3.4 https_proxy=http://1.2.3.4;
fi

if [ "$COVERAGE" ]; then
    COVERAGE_FNAME="/tmp/test_coverage.xml"
    COVERAGE="-s --cov=pandas --cov-report=xml:$COVERAGE_FNAME"
fi

# If no X server is found, we use xvfb to emulate it
if [[ $(uname) == "Linux" && -z $DISPLAY ]]; then
    export DISPLAY=":0"
    XVFB="xvfb-run "
fi

PYTEST_CMD="${XVFB}pytest -m \"$PATTERN\" -n $PYTEST_WORKERS --dist=loadfile -s --strict --durations=30 --junitxml=test-data.xml $TEST_ARGS $COVERAGE pandas"

if [[ $(uname) != "Linux"  && $(uname) != "Darwin" ]]; then
    # GH#37455 windows py38 build appears to be running out of memory
    #  run tests in multiple phases to limit memory footprint
    PYTEST_CMD="$PYTEST_CMD/tests/window"
    PYTEST_CMD="$PYTEST_CMD/tests/api"
    PYTEST_CMD="$PYTEST_CMD/tests/arithmetic"
    PYTEST_CMD="$PYTEST_CMD/tests/arrays"
    PYTEST_CMD="$PYTEST_CMD/tests/base"
    PYTEST_CMD="$PYTEST_CMD/tests/computation"
    PYTEST_CMD="$PYTEST_CMD/tests/config"
    PYTEST_CMD="$PYTEST_CMD/tests/dtypes"
    PYTEST_CMD="$PYTEST_CMD/tests/extension"
    PYTEST_CMD="$PYTEST_CMD/tests/frame"
    PYTEST_CMD="$PYTEST_CMD/tests/generic"
    PYTEST_CMD="$PYTEST_CMD/tests/groupby"
    PYTEST_CMD="$PYTEST_CMD/tests/indexes"
    PYTEST_CMD="$PYTEST_CMD/tests/indexing"
    PYTEST_CMD="$PYTEST_CMD/tests/internals"
    PYTEST_CMD="$PYTEST_CMD/tests/io"
    PYTEST_CMD="$PYTEST_CMD/tests/libs"
    PYTEST_CMD="$PYTEST_CMD/tests/plotting"
    PYTEST_CMD="$PYTEST_CMD/tests/reductions"
    PYTEST_CMD="$PYTEST_CMD/tests/resample"
    PYTEST_CMD="$PYTEST_CMD/tests/reshape"
    PYTEST_CMD="$PYTEST_CMD/tests/scalar"
    PYTEST_CMD="$PYTEST_CMD/tests/series"
    PYTEST_CMD="$PYTEST_CMD/tests/tools"
    PYTEST_CMD="$PYTEST_CMD/tests/tseries"
    PYTEST_CMD="$PYTEST_CMD/tests/tslibs"
    PYTEST_CMD="$PYTEST_CMD/tests/util"
    PYTEST_CMD="$PYTEST_CMD/tests/window"
    PYTEST_CMD="$PYTEST_CMD/tests/test*.py"
fi

echo $PYTEST_CMD
sh -c "$PYTEST_CMD"

if [[ "$COVERAGE" && $? == 0 && "$TRAVIS_BRANCH" == "master" ]]; then
    echo "uploading coverage"
    echo "bash <(curl -s https://codecov.io/bash) -Z -c -f $COVERAGE_FNAME"
          bash <(curl -s https://codecov.io/bash) -Z -c -f $COVERAGE_FNAME
fi
