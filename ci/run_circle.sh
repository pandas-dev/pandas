#!/usr/bin/env bash

echo "[running tests]"
export PATH="$MINICONDA_DIR/bin:$PATH"

source activate pandas

echo "pytest --strict --junitxml=$CIRCLE_TEST_REPORTS/reports/junit.xml $@ pandas"
pytest -W error:ResourceWarning --strict --junitxml=$CIRCLE_TEST_REPORTS/reports/junit.xml $@ pandas
