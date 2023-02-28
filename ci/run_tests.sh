#!/bin/bash -e

# Workaround for pytest-xdist (it collects different tests in the workers if PYTHONHASHSEED is not set)
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
export PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')

# May help reproduce flaky CI builds if set in subsequent runs
echo PYTHONHASHSEED=$PYTHONHASHSEED

if [[ "not network" == *"$PATTERN"* ]]; then
    export http_proxy=http://1.2.3.4 https_proxy=http://1.2.3.4;
fi

if [[ "$COVERAGE" == "true" ]]; then
    COVERAGE="-s --cov=pandas --cov-report=xml --cov-append"
else
    COVERAGE="" # We need to reset this for COVERAGE="false" case
fi

# If no X server is found, we use xvfb to emulate it
if [[ $(uname) == "Linux" && -z $DISPLAY ]]; then
    export DISPLAY=":0"
    XVFB="xvfb-run "
fi

PYTEST_CMD="${XVFB}pytest -r fEs -n $PYTEST_WORKERS --dist=loadfile $TEST_ARGS $COVERAGE $PYTEST_TARGET"

if [[ "$PATTERN" ]]; then
  PYTEST_CMD="$PYTEST_CMD -m \"$PATTERN\""
fi

echo $PYTEST_CMD
sh -c "$PYTEST_CMD"

if [[ "$PANDAS_DATA_MANAGER" != "array" && "$PYTEST_TARGET" == "pandas" ]]; then
    # The ArrayManager tests should have already been run by PYTEST_CMD if PANDAS_DATA_MANAGER was already set to array
    # If we're targeting specific files, e.g. test_downstream.py, don't run.
    PYTEST_AM_CMD="PANDAS_DATA_MANAGER=array pytest -n $PYTEST_WORKERS --dist=loadfile $TEST_ARGS $COVERAGE pandas"

    if [[ "$PATTERN" ]]; then
      PYTEST_AM_CMD="$PYTEST_AM_CMD -m \"$PATTERN and arraymanager\""
    else
      PYTEST_AM_CMD="$PYTEST_AM_CMD -m \"arraymanager\""
    fi

    echo $PYTEST_AM_CMD
    sh -c "$PYTEST_AM_CMD"
fi
