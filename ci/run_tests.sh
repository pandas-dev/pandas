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

# Travis does not have have an X server
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    DISPLAY=DISPLAY=:99.0
    if [ `uname -m` = 'aarch64' ]; then
        PYTEST_CMD="xvfb-run -e /dev/stdout pytest -m \"$PATTERN\" -s --strict --durations=10 --junitxml=test-data.xml $TEST_ARGS $COVERAGE pandas"
    else
        PYTEST_CMD="xvfb-run -e /dev/stdout $PYTEST_CMD"
    fi
fi

PYTEST_CMD="${XVFB}pytest -m \"$PATTERN\" -n auto --dist=loadfile -s --strict --durations=10 --junitxml=test-data.xml $TEST_ARGS $COVERAGE pandas"

echo $PYTEST_CMD
if [ `uname -m` = 'aarch64' ]; then
   sudo sh -c "$PYTEST_CMD"
else
    sh -c "$PYTEST_CMD"
fi

if [[ "$COVERAGE" && $? == 0 && "$TRAVIS_BRANCH" == "master" ]]; then
    echo "uploading coverage"
    echo "bash <(curl -s https://codecov.io/bash) -Z -c -F $TYPE -f $COVERAGE_FNAME"
          bash <(curl -s https://codecov.io/bash) -Z -c -f $COVERAGE_FNAME
fi
