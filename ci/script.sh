#!/bin/bash

echo "inside $0"

python setup.py build_ext install

# HACK: pandas is a statsmodels dependency
# and pandas has a reverse dep (for tests only)
# so we need to install it after pandas
if [ x"$FULL_DEPS" == x"true" ]; then
    # we need at least 0.5.0-dev, just pick a commit for now
    pip install git+git://github.com/statsmodels/statsmodels@c9062e43b8a5f7385537ca95#egg=statsmod

fi;

if [ x"$VBENCH" != x"true" ]; then
    nosetests --exe -w /tmp -A "not slow" pandas;
    exit
fi
if [ x"$VBENCH" == x"true" ]; then
    python vb_suite/perf_HEAD.py;
    exit
fi
