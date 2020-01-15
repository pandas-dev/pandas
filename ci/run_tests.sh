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

if [[ "not clipboard" != *"$PATTERN"* ]]; then
    echo "DISPLAY (original): $DISPLAY"

    # An X server has to exist, and DISPLAY set, for the clipboard (and its tests) to work
    if [ -z $DISPLAY ]; then
        export DISPLAY=":0"
        XVFB="xvfb-run -e /dev/stdout "
    fi

    echo "xsel path: $(which xsel)"
    echo "DISPLAY: $DISPLAY"

    echo "testing xsel from terminal"
    xvfb-run -e /dev/stdout "echo \"terminal clipboard works\" | xsel -i ; echo \"clipboard content: $(xsel -o)\""

    echo "testing xsel from pandas"
    python -c "import pandas.io.clipboard ; pandas.io.clipboard.clipboard_set('pandas clipboard works') ; print(f'clipboard content: {pandas.io.clipboard.clipboard_get()}')"
fi

PYTEST_CMD="${XVFB}pytest -m \"$PATTERN\" -n auto --dist=loadfile -s --strict --durations=10 --junitxml=test-data.xml $TEST_ARGS $COVERAGE pandas/tests/io/test_clipboard.py"

echo $PYTEST_CMD
sh -c "$PYTEST_CMD"

if [[ "$COVERAGE" && $? == 0 && "$TRAVIS_BRANCH" == "master" ]]; then
    echo "uploading coverage"
    echo "bash <(curl -s https://codecov.io/bash) -Z -c -f $COVERAGE_FNAME"
          bash <(curl -s https://codecov.io/bash) -Z -c -f $COVERAGE_FNAME
fi
