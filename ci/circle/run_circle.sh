#!/usr/bin/env bash

echo "[running tests]"
export PATH="$MINICONDA_DIR/bin:$PATH"

source activate pandas

echo "pytest --strict --durations=10 --color=no --junitxml=$CIRCLE_TEST_REPORTS/reports/junit.xml $@ pandas"
pytest       --strict --durations=10 --color=no --junitxml=$CIRCLE_TEST_REPORTS/reports/junit.xml $@ pandas
