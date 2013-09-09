#!/bin/sh

CURDIR=$(pwd)
BASEDIR=$(cd "$(dirname "$0")"; pwd)
python "$BASEDIR"/vb_suite/test_perf.py $@
