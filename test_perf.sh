#!/bin/sh

CURDIR=$(pwd)
BASEDIR=$(readlink -f $(dirname $0  ))

echo "Use vbench to compare the performance of one commit against another."
echo "Make sure the python 'vbench' library is installed..\n"

cd "$BASEDIR/vb_suite/"
python test_perf.py $@

cd "$CURDIR"
