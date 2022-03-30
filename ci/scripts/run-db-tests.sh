#!/bin/bash

set -e

python setup.py build_ext --inplace -j4

/etc/init.d/postgresql start
python -m pytest pandas/tests/io/test_sql.py
