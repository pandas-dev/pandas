#!/bin/bash -e

# Workaround for pytest-xdist (it collects different tests in the workers if PYTHONHASHSEED is not set)
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')
export PYTHONHASHSEED

COVERAGE="-s --cov=pandas --cov-report=xml --cov-append --cov-config=pyproject.toml"

PYTEST_CMD="MESONPY_EDITABLE_VERBOSE=1 PYTHONDEVMODE=1 PYTHONWARNDEFAULTENCODING=1 pytest -r fE -n $PYTEST_WORKERS --dist=worksteal $TEST_ARGS $COVERAGE $PYTEST_TARGET"

if [[ "$PATTERN" ]]; then
  PYTEST_CMD="$PYTEST_CMD -m \"$PATTERN\""
fi

echo "$PYTEST_CMD"
sh -c "$PYTEST_CMD"
