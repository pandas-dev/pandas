#!/bin/bash -e

# Workaround for pytest-xdist (it collects different tests in the workers if PYTHONHASHSEED is not set)
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
export PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')

# May help reproduce flaky CI builds if set in subsequent runs
echo PYTHONHASHSEED=$PYTHONHASHSEED

COVERAGE="-s --cov=pandas --cov-report=xml --cov-append --cov-config=pyproject.toml"

PYTEST_CMD="MESONPY_EDITABLE_VERBOSE=1 PYTHONDEVMODE=1 PYTHONWARNDEFAULTENCODING=1 pytest -r fE -n $PYTEST_WORKERS --dist=worksteal $TEST_ARGS $COVERAGE $PYTEST_TARGET"

if [[ "$PATTERN" ]]; then
  PYTEST_CMD="$PYTEST_CMD -m \"$PATTERN and not bodo_udf_engine\""
else
  PYTEST_CMD="$PYTEST_CMD -m \"not bodo_udf_engine\""
fi

echo $PYTEST_CMD
sh -c "$PYTEST_CMD"

# Workaround for running bodo tests. Needs to be in a separate session to prevent
# conflicts with numba extensions and run without PYTHONDEVMODE=1 since it can cause segmentation faults during compilation.
if [[ "$PYTEST_WORKERS" == "0" ]]; then
  PYTEST_CMD_BODO_UDF_ENGINE="MESONPY_EDITABLE_VERBOSE=1 PYTHONWARNDEFAULTENCODING=1 pytest -r fE -n $PYTEST_WORKERS --dist=worksteal $TEST_ARGS $COVERAGE $PYTEST_TARGET -m \"bodo_udf_engine\""
  echo "Running Bodo Tests..."
  echo $PYTEST_CMD_BODO_UDF_ENGINE
  sh -c "$PYTEST_CMD_BODO_UDF_ENGINE"
fi
