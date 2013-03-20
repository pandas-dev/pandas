#!/bin/sh

CURDIR=$(pwd)
BASEDIR=$(readlink -f $(dirname $0  ))

echo "This script compares the performance of two commits."
echo "Make sure the python 'vbench' library is installed.\n"
echo "Setting the BUILD_CACHE_DIR env var to a temp directory will"
echo "potentially speed up subsequent runs.\n"


cd "$BASEDIR/vb_suite/"
python test_perf.py $@

cd "$CURDIR"
