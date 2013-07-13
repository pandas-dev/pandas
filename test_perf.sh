#!/bin/sh

CURDIR=$(pwd)
BASEDIR=$(readlink -f $(dirname $0  ))

python "$BASEDIR"/vb_suite/test_perf.py $@
