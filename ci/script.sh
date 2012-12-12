#!/bin/bash

echo "inside $0"

if [ x"$VBENCH" != x"true" ]; then
    nosetests --exe -w /tmp -A "not slow" pandas;
    exit
fi
if [ x"$VBENCH" == x"true" ]; then
    python vb_suite/perf_HEAD.py;
    exit
fi
