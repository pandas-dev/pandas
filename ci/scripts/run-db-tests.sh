#!/bin/bash

set -e

python setup.py build_ext --inplace -j2
python -m pytest pandas/tests/io/test_sql.py
