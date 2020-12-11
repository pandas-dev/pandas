#!/bin/bash -e
# Special build for python3.9 until numpy puts its own wheels up

pip install --no-deps -U pip wheel setuptools
pip install cython numpy python-dateutil pytz pytest pytest-xdist hypothesis

python setup.py build_ext -inplace
python -m pip install --no-build-isolation -e .

python -c "import sys; print(sys.version_info)"
python -c "import pandas as pd"
python -c "import hypothesis"
