#!/bin/bash

if [ -z "$COVERAGE" ]; then
   echo "coverage is not selected for this build"
   exit 0
fi

source activate pandas

echo "uploading coverage"
bash <(curl -s https://codecov.io/bash) -Z -c -F single -f /tmp/cov-single.xml
bash <(curl -s https://codecov.io/bash) -Z -c -F multiple -f /tmp/cov-multiple.xml
