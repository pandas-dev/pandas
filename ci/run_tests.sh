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

echo $PYTEST_CMD

if [[ $(uname) != "Linux"  && $(uname) != "Darwin" ]]; then
    # GH#37455 windows py38 build appears to be running out of memory
    #  run tests in multiple phases to limit memory footprint
    sh -c "$PYTEST_CMD/tests/a*"
    sh -c "$PYTEST_CMD/tests/b*"
    sh -c "$PYTEST_CMD/tests/c*"
    sh -c "$PYTEST_CMD/tests/d*"
    sh -c "$PYTEST_CMD/tests/e*"
    sh -c "$PYTEST_CMD/tests/g*"
    sh -c "$PYTEST_CMD/tests/i*"
    sh -c "$PYTEST_CMD/tests/l*"
    sh -c "$PYTEST_CMD/tests/p*"
    sh -c "$PYTEST_CMD/tests/r*"
    sh -c "$PYTEST_CMD/tests/s*"
    sh -c "$PYTEST_CMD/tests/t*"
    sh -c "$PYTEST_CMD/tests/u*"
    sh -c "$PYTEST_CMD/tests/w*"
else
    sh -c "$PYTEST_CMD"
fi

echo $PYTEST_CMD
sh -c "$PYTEST_CMD"

if [[ "$COVERAGE" && $? == 0 && "$TRAVIS_BRANCH" == "master" ]]; then
    echo "uploading coverage"
    echo "bash <(curl -s https://codecov.io/bash) -Z -c -f $COVERAGE_FNAME"
          bash <(curl -s https://codecov.io/bash) -Z -c -f $COVERAGE_FNAME
fi
