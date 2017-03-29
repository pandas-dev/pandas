#!/bin/bash

if [ -z "$COVERAGE" ]; then
   echo "no upload of coverage is needed"
   exit 0
fi

source activate pandas

codecov --file -c -F single /tmp/cov-single.xml
codecov --file -c -F multiple /tmp/cov-multiple.xml
