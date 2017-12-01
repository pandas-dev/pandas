#!/bin/bash

# Workaround for pytest-xdist flaky collection order
# https://github.com/pytest-dev/pytest/issues/920
# https://github.com/pytest-dev/pytest/issues/1075
export PYTHONHASHSEED=$(python -c 'import random; print(random.randint(1, 4294967295))')

pytest pandas --skip-slow --skip-network -m "not single" -n 4 -r sxX --strict "$@"
